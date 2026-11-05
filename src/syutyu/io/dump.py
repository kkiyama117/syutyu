"""LAMMPS の custom dump を MDAnalysis の Universe として読む。
"""
from pathlib import Path

import MDAnalysis
import numpy as np

COORDINATE_COLUMNS = ((('xu', 'yu', 'zu'), 'unwrapped'), (('x', 'y', 'z'), 'unscaled'),
                      (('xs', 'ys', 'zs'), 'scaled'))


def _coordinate_columns(columns):
    """座標の列と MDAnalysis に渡す規約を選ぶ。"""
    for names, convention in COORDINATE_COLUMNS:
        if all(name in columns for name in names):
            return names, convention
    raise ValueError('Need x/y/z, xs/ys/zs, or xu/yu/zu coordinates')


def _scan(path):
    """フレームを先頭から走査し、列の並びとフレームごとのセル下限を返す。

    数値に直すのは step・原子数・セル境界だけで、原子行は列数を見るだけにしてある。
    最初の違反で止まるので、MDAnalysis を建てた後の検査より先にこのメッセージが出る。
    """
    columns, origins = None, []
    count0, last_step = None, None
    with Path(path).open(encoding='utf-8') as stream:
        while first := stream.readline():
            if first.strip() != 'ITEM: TIMESTEP':
                raise ValueError('Expected ITEM: TIMESTEP')
            step = int(stream.readline())
            if last_step is not None and step <= last_step:
                raise ValueError('Need increasing timesteps; split trajectories at reset_timestep')
            last_step = step
            if stream.readline().strip() != 'ITEM: NUMBER OF ATOMS':
                raise ValueError('Missing atom-count header')
            count = int(stream.readline())
            if count < 1:
                raise ValueError('Empty atom set')
            if count0 is None:
                count0 = count
            elif count != count0:
                raise ValueError('Atom count changed across frames')
            if stream.readline().strip() != 'ITEM: BOX BOUNDS pp pp pp':
                raise ValueError('Only orthogonal periodic pp pp pp cells are supported (no triclinic tilt)')
            bounds = np.array([list(map(float, stream.readline().split())) for _ in range(3)])
            if bounds.shape != (3, 2) or not np.isfinite(bounds).all():
                raise ValueError('Invalid box bounds')
            if np.any(bounds[:, 1] - bounds[:, 0] <= 0):
                raise ValueError('Box lengths must be positive')
            header = stream.readline().split()
            if header[:2] != ['ITEM:', 'ATOMS']:
                raise ValueError('Missing atom-column header')
            frame_columns = header[2:]
            if columns is None:
                columns = frame_columns
            if (frame_columns != columns or len(set(frame_columns)) != len(frame_columns)
                    or 'id' not in frame_columns):
                raise ValueError('Need stable, unique atom columns including id')
            present_velocity = [name in frame_columns for name in ('vx', 'vy', 'vz')]
            if any(present_velocity) and not all(present_velocity):
                raise ValueError('Velocity output must include all of vx/vy/vz or none of them')
            rows = [stream.readline().split() for _ in range(count)]
            if any(len(row) != len(frame_columns) for row in rows):
                raise ValueError('Incomplete atom records')
            _coordinate_columns(frame_columns)
            origins.append(bounds[:, 0])
    if not origins:
        raise ValueError('Need increasing timesteps; split trajectories at reset_timestep')
    return columns, np.array(origins)


def _shift(offsets):
    """MDAnalysis が引いたセル下限を足し戻す変換を作る。"""
    def transform(timestep):
        timestep.positions = timestep.positions + offsets[timestep.frame]
        return timestep
    return transform


def read_trajectory(path):
    """LAMMPS の custom dump を読み、MDAnalysis の Universe を返す。

    Args:
        path (str | os.PathLike): custom dump のファイル。``ITEM: TIMESTEP`` で始まるフレームを並べたもの。

    Returns:
        MDAnalysis.Universe: 座標は LAMMPS の絶対座標、セルは辺の長さ。``mol`` 列がなければ resid は全て 0、``type`` と ``element`` 列がなければ空文字になる。速度は vx/vy/vz があれば入る。

    Raises:
        ValueError: 直方体・3方向周期（``pp pp pp``）以外のセル、重複 ID、原子数がフレーム間で変わる、step の減少、非有限値、vx/vy/vz の一部だけの出力。

    Note:
        座標と速度は float32、座標は LAMMPS の絶対座標。xu/yu/zu か、x/y/z と image flag があれば非折り返し座標になり、xs/ys/zs は image flag があっても折り返したまま。``q`` などの追加列の扱いは MDAnalysis に従う（``q`` は ``atoms.charges``）。単位は入力のまま扱い、``reset_timestep`` をまたぐ軌跡は分割して読む。
    """
    columns, origins = _scan(path)
    _names, convention = _coordinate_columns(columns)
    unwrap_images = convention == 'unscaled' and all(name in columns for name in ('ix', 'iy', 'iz'))
    universe = MDAnalysis.Universe(str(path), format='LAMMPSDUMP', dt=1.0,
                                   lammps_coordinate_convention=convention,
                                   unwrap_images=unwrap_images)
    universe.trajectory.add_transformations(_shift(origins if convention != 'scaled' else origins * 2))
    ids = universe.atoms.ids
    if len(set(ids.tolist())) != len(ids):
        raise ValueError('Duplicate atom IDs')
    for timestep in universe.trajectory:
        if not np.isfinite(timestep.positions).all():
            raise ValueError('Non-finite coordinates')
        if timestep.has_velocities and not np.isfinite(timestep.velocities).all():
            raise ValueError('Non-finite velocities')
    n_atoms = universe.atoms.n_atoms
    if 'mol' not in columns:
        universe.residues.resids = [0]
    if 'type' not in columns:
        universe.add_TopologyAttr('type', np.full(n_atoms, '', dtype=object))
    if 'element' not in columns:
        universe.add_TopologyAttr('element', np.full(n_atoms, '', dtype=object))
    universe.trajectory[0]
    return universe


def _steps(universe):
    """1フレームに1つ、LAMMPS の step を並べた整数配列を返す。"""
    steps = [int(timestep.data['step']) if 'step' in timestep.data else int(round(timestep.time))
             for timestep in universe.trajectory]
    return np.array(steps, dtype=np.int64)


def _frames(universe, *, start_step=None, stop_step=None, stride=1):
    """step の閉区間に入るフレーム位置を、間引きながら返す。"""
    if not isinstance(stride, int) or isinstance(stride, bool) or stride < 1:
        raise ValueError('Frame stride must be a positive integer')
    steps = _steps(universe)
    lower = steps[0] if start_step is None else start_step
    upper = steps[-1] if stop_step is None else stop_step
    if lower > upper:
        raise ValueError('start_step must not exceed stop_step')
    selected = np.flatnonzero((steps >= lower) & (steps <= upper))[::stride]
    if not len(selected):
        raise ValueError('No trajectory frames selected')
    return selected


def _select(universe, expression):
    """MDAnalysis の選択式を、原子数と同じ長さの真偽マスクにする。"""
    try:
        group = universe.select_atoms(expression)
    except (MDAnalysis.exceptions.SelectionError, MDAnalysis.exceptions.NoDataError,
            TypeError, AttributeError) as error:
        raise ValueError(f'Invalid selection {expression!r}: {error}') from None
    if not len(group):
        raise ValueError(f'Selection matches no atoms: {expression}')
    mask = np.zeros(universe.atoms.n_atoms, dtype=bool)
    mask[group.indices] = True
    return mask


def _positions(universe, frames):
    """指定フレームの座標を frames × atoms × 3 の float64 配列で返す。"""
    return np.array([np.asarray(universe.trajectory[int(frame)].positions, dtype=float)
                     for frame in frames])


def _boxes(universe, frames):
    """指定フレームのセルの辺の長さを frames × 3 の float64 配列で返す。"""
    boxes = []
    for frame in frames:
        dimensions = universe.trajectory[int(frame)].dimensions
        if dimensions is None:
            raise ValueError('Universe has no cell dimensions')
        if not np.allclose(dimensions[3:], 90, rtol=0, atol=1e-6):
            raise ValueError('Only orthogonal cells are supported')
        boxes.append(np.asarray(dimensions[:3], dtype=float))
    return np.array(boxes)


def _velocities(universe, frames):
    """指定フレームの速度を frames × atoms × 3 の float64 配列で返す。"""
    if not universe.trajectory.ts.has_velocities:
        return None
    return np.array([np.asarray(universe.trajectory[int(frame)].velocities, dtype=float)
                     for frame in frames])


def _fixed_box(boxes):
    """全フレームでセルの辺の長さが変わっていないかを返す。"""
    return np.allclose(boxes, boxes[0], rtol=0, atol=1e-10)


def _source_dump(universe):
    """結果に残す読み込み元のパスを返す。"""
    filename = getattr(universe, 'filename', None)
    return None if filename is None else str(Path(filename).resolve())
