"""力場に依存しない原子 MSD。
"""
from collections.abc import Mapping

import numpy as np

from ..io.dump import frame_indices
from ..io.store import Result
from ..kernels import msd_multi_origin, unwrap


def _mass_per_atom(trajectory, com_masses, laboratory):
    """`com_masses` を検証し、原子ごとの質量 (g/mol) と使った type ごとの質量を返す。

    引数:
        - `trajectory` — `Trajectory`。`type` ラベルを読む。
        - `com_masses` — 検証する値。`None` では呼ばない。
        - `laboratory` — 呼び出し側の指定。同時指定を拒否するために見る。
    戻り値:
        原子ごとの質量の配列と、使われている type の質量の辞書（キーは昇順）。
    例外:
        `ValueError` — マッピングでない、`laboratory=True` と同時、type 列がない、質量がない、質量が有限の正の数でない。
    """
    if not isinstance(com_masses, Mapping):
        raise ValueError('com_masses must be a mapping from type label to mass')
    if laboratory:
        raise ValueError('com_masses cannot be combined with laboratory=True')
    if 'type' not in trajectory.labels:
        raise ValueError('com_masses needs a type column in the dump')
    present = sorted(set(trajectory.labels['type'].tolist()))
    missing = [label for label in present if label not in com_masses]
    if missing:
        raise ValueError(f'com_masses has no mass for types: {missing}')
    mass_of = {}
    for label in present:
        try:
            mass = float(com_masses[label])
        except (TypeError, ValueError):
            raise ValueError(f'com_masses: mass of type {label} must be finite and > 0') from None
        if not np.isfinite(mass) or mass <= 0:
            raise ValueError(f'com_masses: mass of type {label} must be finite and > 0')
        mass_of[label] = mass
    return np.array([mass_of[label] for label in trajectory.labels['type']]), mass_of


def mean_squared_displacement(trajectory, selections, *, dt_ps, start_step=None, stop_step=None,
                              stride=1, laboratory=False, origin_mode='multiple', com_masses=None):
    """選択ごとに原子 MSD を計算する。拡散フィットは行わない。

    引数:
        - `trajectory` — `Trajectory`。
        - `selections` — 選択式の並び。各選択が1つの系列になる。
        - `dt_ps` — LAMMPS の 1 step 当たりの ps。
        - `start_step`、`stop_step`、`stride` — 使うフレームの範囲と間引き。等間隔であること。
        - `laboratory` — 可変セルで、セル伸縮を含む実験室座標の MSD を明示的に求める。
        - `origin_mode` — `multiple`（全開始時刻）か `single`（最初の選択フレーム）。
        - `com_masses` — dump の type ラベル → 質量 g/mol。与えると全原子の質量重み付き重心の変位を各原子から引く。選択ごとの重心は引かない。
    戻り値:
        `Result`。`tables['msd']` に `lag_ps`、`MSD_<i>_A2`、`origin_atom_pairs_<i>`。`summary` に `com_drift_removed`。除去したときは `com_masses` と `com_max_displacement_A`（Å）も入る。
    例外:
        `ValueError` — `dt_ps` が正でない、原点モードが不正、フレームが2未満か等間隔でない、選択が空、可変セルなのに非折り返し座標と `laboratory=True` がそろわない。
        `ValueError` — `com_masses` がマッピングでない、type 列がない、使われている type の質量がない、質量が有限の正の数でない、`laboratory=True` と同時に与えた。
    注意:
        原子の MSD であり、分子重心の MSD ではない。重なった時間原点は独立でなく、平均二乗変位から平衡を主張しない。
    """
    if not np.isfinite(dt_ps) or dt_ps <= 0:
        raise ValueError('dt_ps must be positive ps per LAMMPS timestep')
    if origin_mode not in ('multiple', 'single'):
        raise ValueError('origin_mode must be multiple or single')
    frames = frame_indices(trajectory, start_step=start_step, stop_step=stop_step, stride=stride)
    steps = trajectory.steps[frames]
    if len(frames) < 2 or not np.all(np.diff(steps) == steps[1] - steps[0]):
        raise ValueError('MSD needs at least two uniformly spaced selected frames')
    if not selections:
        raise ValueError('Specify at least one atom selection')
    if not trajectory.fixed_box and (not laboratory or trajectory.unwrapped is None):
        raise ValueError('Variable-cell MSD requires xu/yu/zu or image flags AND laboratory=True; '
                         'wrapped NPT reconstruction is intentionally unsupported')
    masses = mass_of = None
    if com_masses is not None:
        masses, mass_of = _mass_per_atom(trajectory, com_masses, laboratory)
    if trajectory.unwrapped is not None:
        coordinates = trajectory.unwrapped[frames]
        source = trajectory.unwrap_source
    else:
        # Unwrap BEFORE frame thinning; otherwise stride can hide boundary crossings.
        coordinates = unwrap(trajectory.positions, trajectory.boxes[0])[frames]
        source = 'minimum-image increments (conditional)'
    com_shift = None
    if masses is not None:
        com = np.sum(coordinates * masses[None, :, None], axis=1) / masses.sum()
        com_shift = com - com[0]
        coordinates = coordinates - com_shift[:, None, :]
    table = {'lag_ps': np.arange(len(frames)) * (steps[1] - steps[0]) * dt_ps}
    info = []
    for index, expression in enumerate(selections):
        mask = trajectory.select(expression)
        if origin_mode == 'multiple':
            msd, samples = msd_multi_origin(coordinates[:, mask])
        else:
            displacement = coordinates[:, mask] - coordinates[:1, mask]
            msd = np.mean(np.sum(displacement ** 2, axis=-1), axis=1)
            samples = np.full(len(frames), int(mask.sum()), dtype=int)
        key = f'MSD_{index}_A2'
        table[key] = msd
        table[f'origin_atom_pairs_{index}'] = samples
        info.append({'selection': expression, 'column': key, 'atoms': int(mask.sum())})
    summary = {
        'kind': 'msd', 'title': 'Species/site atomic MSD', 'atoms': len(trajectory.ids),
        'method': 'custom', 'origin_mode': origin_mode, 'origin_step': int(steps[0]),
        'source_dump': trajectory.source_path,
        'selected_frames': len(frames), 'selected_steps': [int(steps[0]), int(steps[-1])],
        'dt_ps': dt_ps, 'frame_stride': stride, 'msd_selections': info,
        'coordinate_frame': 'laboratory (includes affine box motion)' if not trajectory.fixed_box else 'fixed box',
        'unwrap_source': source,
        'com_drift_removed': masses is not None,
    }
    limitations = ['Atomic MSD, not molecular center-of-mass MSD; no COM drift removal',
                   'Overlapping time origins are correlated; no diffusion coefficient inferred',
                   'Equilibrium not confirmed; variable-box laboratory MSD is not a corrected diffusion estimate']
    if masses is not None:
        summary['com_masses'] = mass_of
        summary['com_max_displacement_A'] = float(np.max(np.linalg.norm(com_shift, axis=1)))
        limitations[0] = ('Atomic MSD, not molecular center-of-mass MSD; '
                          'mass-weighted system COM drift removed')
    summary['limitations'] = limitations
    return Result({'msd': table}, summary)
