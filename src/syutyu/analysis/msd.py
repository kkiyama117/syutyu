"""力場に依存しない原子 MSD。
"""
import numpy as np

from ..io.dump import frame_indices
from ..io.store import Result
from ..kernels import msd_multi_origin, unwrap


def mean_squared_displacement(trajectory, selections, *, dt_ps, start_step=None, stop_step=None,
                              stride=1, laboratory=False, origin_mode='multiple'):
    """選択ごとに原子 MSD を計算する。拡散フィットや重心除去は行わない。

    引数:
        - `trajectory` — `Trajectory`。
        - `selections` — 選択式の並び。各選択が1つの系列になる。
        - `dt_ps` — LAMMPS の 1 step 当たりの ps。
        - `start_step`、`stop_step`、`stride` — 使うフレームの範囲と間引き。等間隔であること。
        - `laboratory` — 可変セルで、セル伸縮を含む実験室座標の MSD を明示的に求める。
        - `origin_mode` — `multiple`（全開始時刻）か `single`（最初の選択フレーム）。
    戻り値:
        `Result`。`tables['msd']` に `lag_ps`、`MSD_<i>_A2`、`origin_atom_pairs_<i>`。
    例外:
        `ValueError` — `dt_ps` が正でない、原点モードが不正、フレームが2未満か等間隔でない、選択が空、可変セルなのに非折り返し座標と `laboratory=True` がそろわない。
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
    if trajectory.unwrapped is not None:
        coordinates = trajectory.unwrapped[frames]
        source = trajectory.unwrap_source
    else:
        # Unwrap BEFORE frame thinning; otherwise stride can hide boundary crossings.
        coordinates = unwrap(trajectory.positions, trajectory.boxes[0])[frames]
        source = 'minimum-image increments (conditional)'
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
    return Result({'msd': table}, {
        'kind': 'msd', 'title': 'Species/site atomic MSD', 'atoms': len(trajectory.ids),
        'method': 'custom', 'origin_mode': origin_mode, 'origin_step': int(steps[0]),
        'source_dump': trajectory.source_path,
        'selected_frames': len(frames), 'selected_steps': [int(steps[0]), int(steps[-1])],
        'dt_ps': dt_ps, 'frame_stride': stride, 'msd_selections': info,
        'coordinate_frame': 'laboratory (includes affine box motion)' if not trajectory.fixed_box else 'fixed box',
        'unwrap_source': source,
        'limitations': ['Atomic MSD, not molecular center-of-mass MSD; no COM drift removal',
                        'Overlapping time origins are correlated; no diffusion coefficient inferred',
                        'Equilibrium not confirmed; variable-box laboratory MSD is not a corrected diffusion estimate'],
    })
