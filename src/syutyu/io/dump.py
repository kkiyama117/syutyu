"""LAMMPS の custom dump を読み、種別・元素・分子で原子を明示的に選ぶ。
"""
from dataclasses import dataclass
from pathlib import Path
import numpy as np


@dataclass
class Trajectory:
    """dump 1つ分のフレーム列と、原子ごとの固定ラベル。

    Attributes:
        steps (numpy.ndarray): 各フレームの LAMMPS step（整数）。
        ids (numpy.ndarray): 原子 ID。全フレームで同じ昇順にそろえてある。
        positions (numpy.ndarray): 各フレームの座標。Å、周期セル内の絶対座標。
        boxes (numpy.ndarray): 各フレームの直方体セルの辺の長さ。Å。
        origins (numpy.ndarray): 各フレームのセル原点。Å。
        labels (dict): dump にあった ``type``・``element``・``mol`` の原子ごとの値（文字列）。
        unwrapped (numpy.ndarray | None): 明示的な非折り返し座標。dump に xu/yu/zu か image flag があるときだけ入り、なければ None。
        unwrap_source (str | None): ``unwrapped`` の作り方（``xu/yu/zu`` か ``image flags``）。なければ None。
        source_path (str | None): 読み込んだファイルの絶対パス。
        velocities (numpy.ndarray | None): 各フレームの速度。ID 順、保存時の単位のまま。vx/vy/vz がそろわなければ None。
    """
    steps: np.ndarray
    ids: np.ndarray
    positions: np.ndarray       # Cartesian Å, absolute box coordinates
    boxes: np.ndarray           # orthogonal lengths per frame
    origins: np.ndarray
    labels: dict               # static per-atom type, element, mol (strings)
    unwrapped: np.ndarray | None
    unwrap_source: str | None
    source_path: str | None = None
    velocities: np.ndarray | None = None   # frames × atoms × 3, sorted by id; unit as saved

    def select(self, expression):
        """選択式に合う原子の真偽マスクを返す。

        Args:
            expression (str): ``all``、``type:3``、``element:Na``、``mol:12`` のいずれか。ラベルは dump にあった列だけを使う。

        Returns:
            numpy.ndarray: 原子数と同じ長さの真偽配列。

        Raises:
            ValueError: 形式がキー:値 でない、dump にないキーを指す、または1個も一致しない。

        Note:
            元素から原子タイプや力場部位を推測しない。``mol:0`` は未設定を表すことがあり、選んでも分子として意味があるとは限らない。
        """
        if expression == 'all':
            return np.ones(len(self.ids), dtype=bool)
        key, separator, value = expression.partition(':')
        if not separator or key not in self.labels:
            raise ValueError(f'Unavailable selection {expression!r}; dump fields: {list(self.labels)}')
        mask = self.labels[key] == value
        if not np.any(mask):
            raise ValueError(f'Selection matches no atoms: {expression}')
        return mask

    @property
    def fixed_box(self):
        """全フレームでセルの辺と原点が変わっていないかを返す。

        Returns:
            bool: 変化が許容差 ``1e-10`` 未満なら True。可変体積の判定に使う。
        """
        return (np.allclose(self.boxes, self.boxes[0], rtol=0, atol=1e-10)
                and np.allclose(self.origins, self.origins[0], rtol=0, atol=1e-10))


def frame_indices(trajectory, *, start_step=None, stop_step=None, stride=1):
    """step の閉区間に入るフレーム位置を、間引きながら返す。

    Args:
        trajectory (Trajectory): 対象の軌跡。
        start_step (int | None): ``start_step``・``stop_step`` の省略時は最初と最後の step を使う。
        stop_step (int | None): ``start_step``・``stop_step`` の省略時は最初と最後の step を使う。
        stride (int): 選んだフレームを何個おきに残すか（1以上）。

    Returns:
        numpy.ndarray: フレーム位置の整数配列。

    Raises:
        ValueError: ``stride`` が正の整数でない、``start_step > stop_step``、該当フレームがない。

    Note:
        間引きは区間で絞った後に行う。step が等間隔でない軌跡は、MSD や VACF 側の等間隔検査で止まる。
    """
    if not isinstance(stride, int) or isinstance(stride, bool) or stride < 1:
        raise ValueError('Frame stride must be a positive integer')
    lower = trajectory.steps[0] if start_step is None else start_step
    upper = trajectory.steps[-1] if stop_step is None else stop_step
    if lower > upper:
        raise ValueError('start_step must not exceed stop_step')
    selected = np.flatnonzero((trajectory.steps >= lower) & (trajectory.steps <= upper))[::stride]
    if not len(selected):
        raise ValueError('No trajectory frames selected')
    return selected


def read_trajectory(path):
    """LAMMPS の custom dump を読み、Trajectory を返す。

    Args:
        path (str | os.PathLike): custom dump のファイル。``ITEM: TIMESTEP`` で始まるフレームを並べたもの。

    Returns:
        Trajectory: 座標は Å の絶対座標、セルは辺の長さ。速度は vx/vy/vz があれば ID 順で入る。

    Raises:
        ValueError: 直方体・3方向周期（``pp pp pp``）以外のセル、ID やラベルの変化、重複 ID、step の減少、非有限値、vx/vy/vz の一部だけの出力。

    Note:
        xu/yu/zu か image flag があれば非折り返し座標も作り、``unwrap_source`` に残す。``q`` などの追加列は読み捨てる。
        単位は入力のまま扱い、Å 以外への換算はしない。``reset_timestep`` をまたぐ軌跡は分割して読む。
    """
    steps, positions, boxes, origins, unwrapped, velocities = [], [], [], [], [], []
    ids0, labels0, columns0, unwrap_source = None, None, None, None
    with Path(path).open(encoding='utf-8') as stream:
        while first := stream.readline():
            if first.strip() != 'ITEM: TIMESTEP':
                raise ValueError('Expected ITEM: TIMESTEP')
            steps.append(int(stream.readline()))
            if stream.readline().strip() != 'ITEM: NUMBER OF ATOMS':
                raise ValueError('Missing atom-count header')
            n = int(stream.readline())
            if n < 1:
                raise ValueError('Empty atom set')
            if stream.readline().strip() != 'ITEM: BOX BOUNDS pp pp pp':
                raise ValueError('Only orthogonal periodic pp pp pp cells are supported (no triclinic tilt)')
            bounds = np.array([list(map(float, stream.readline().split())) for _ in range(3)])
            if bounds.shape != (3, 2) or not np.isfinite(bounds).all():
                raise ValueError('Invalid box bounds')
            origin, box = bounds[:, 0], bounds[:, 1] - bounds[:, 0]
            if np.any(box <= 0):
                raise ValueError('Box lengths must be positive')
            header = stream.readline().split()
            if header[:2] != ['ITEM:', 'ATOMS']:
                raise ValueError('Missing atom-column header')
            columns = header[2:]
            if columns0 is None:
                columns0 = columns
            if columns != columns0 or len(set(columns)) != len(columns) or 'id' not in columns:
                raise ValueError('Need stable, unique atom columns including id')
            present_velocity = [name in columns for name in ('vx', 'vy', 'vz')]
            if any(present_velocity) and not all(present_velocity):
                raise ValueError('Velocity output must include all of vx/vy/vz or none of them')
            rows = [stream.readline().split() for _ in range(n)]
            if any(len(row) != len(columns) for row in rows):
                raise ValueError('Incomplete atom records')
            data = {name: np.array([row[k] for row in rows]) for k, name in enumerate(columns)}
            ids = data['id'].astype(int)
            order = np.argsort(ids)
            ids = ids[order]
            if len(set(ids)) != n:
                raise ValueError('Duplicate atom IDs')
            labels = {name: data[name][order] for name in ('type', 'element', 'mol') if name in data}
            if ids0 is None:
                ids0, labels0 = ids, labels
            if not np.array_equal(ids, ids0) or any(not np.array_equal(labels[name], labels0[name]) for name in labels):
                raise ValueError('Atom IDs or labels changed across frames')
            def xyz(names):
                return np.column_stack([data[name].astype(float)[order] for name in names])
            explicit = None
            if all(name in data for name in ('xu', 'yu', 'zu')):
                explicit = xyz(('xu', 'yu', 'zu'))
                unwrap_source = 'xu/yu/zu'
            if all(name in data for name in ('x', 'y', 'z')):
                position = xyz(('x', 'y', 'z'))
            elif all(name in data for name in ('xs', 'ys', 'zs')):
                position = origin + xyz(('xs', 'ys', 'zs')) * box
            elif explicit is not None:
                position = explicit.copy()
            else:
                raise ValueError('Need x/y/z, xs/ys/zs, or xu/yu/zu coordinates')
            if explicit is None and all(name in data for name in ('ix', 'iy', 'iz')):
                flags = xyz(('ix', 'iy', 'iz'))
                if not np.array_equal(flags, np.rint(flags)):
                    raise ValueError('Image flags must be integers')
                explicit = position + flags * box
                unwrap_source = 'image flags'
            if not np.isfinite(position).all() or (explicit is not None and not np.isfinite(explicit).all()):
                raise ValueError('Non-finite coordinates')
            if all(present_velocity):
                velocity = xyz(('vx', 'vy', 'vz'))
                if not np.isfinite(velocity).all():
                    raise ValueError('Non-finite velocities')
                velocities.append(velocity)
            positions.append(position)
            boxes.append(box)
            origins.append(origin)
            if explicit is not None:
                unwrapped.append(explicit)
    if not steps or np.any(np.diff(steps) <= 0):
        raise ValueError('Need increasing timesteps; split trajectories at reset_timestep')
    return Trajectory(np.array(steps), ids0, np.array(positions), np.array(boxes), np.array(origins),
                      labels0, np.array(unwrapped) if unwrapped else None, unwrap_source,
                      str(Path(path).resolve()), np.array(velocities) if velocities else None)
