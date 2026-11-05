"""dump の速度から自前で VACF を求める。compute vacf の代替にはしない。
"""
import numpy as np

from ..io.dump import _frames, _select, _source_dump, _steps, _velocities
from ..io.store import Result, _frame


def _correlations(velocities, lags, origin_mode):
    """Average sum_i v_i(m)v_i(m+k) over origins and atoms, per component.

    No cross-atom terms. 'multiple' uses every available origin (Colab);
    'single' fixes origin frame 0 to mirror LAMMPS compute vacf.
    """
    frames, atoms, _ = velocities.shape
    totals = np.zeros((len(lags), 3))
    counts = np.zeros(len(lags), dtype=np.int64)
    if origin_mode == 'multiple':
        for index, lag in enumerate(lags):
            totals[index] = np.sum(velocities[:frames - lag] * velocities[lag:], axis=(0, 1))
            counts[index] = (frames - lag) * atoms
    else:
        for index, lag in enumerate(lags):
            totals[index] = np.sum(velocities[0] * velocities[lag], axis=0)
            counts[index] = atoms
    return totals / counts[:, None], counts


def velocity_autocorrelation(trajectory, selections, *, dt_ps, velocity_unit,
                             start_step=None, stop_step=None, stride=1, origin_mode='multiple',
                             max_lag_ps=None):
    """選択ごとに速度自己相関を計算する。複数時間原点が Colab の定義。

    Args:
        trajectory (MDAnalysis.Universe): 速度を含む軌跡。
        selections (Sequence[str]): MDAnalysis 選択式の並び。各選択が1つの系列になる。
        dt_ps (float): LAMMPS の 1 step 当たりの ps。dump の間隔ではない。
        velocity_unit (str): ``A/fs`` か ``A/ps``。``A/fs`` は二乗する前に換算する。
        start_step (int | None): 使うフレームの範囲と間引き。等間隔であること。
        stop_step (int | None): 使うフレームの範囲と間引き。等間隔であること。
        stride (int): 使うフレームの範囲と間引き。等間隔であること。
        origin_mode (str): ``multiple`` （全開始時刻の平均）か ``single`` （フレーム0のみ）。
        max_lag_ps (float | None): 残す lag の上限（両端を含む）。省略時は全 lag。

    Returns:
        Result: ``tables['vacf']`` に ``lag_ps``、成分 ``VACF_<i>_x_A2_ps2`` など、``VACF_<i>_A2_ps2``、``VACF_<i>_normalized``、``origin_atom_pairs_<i>``。

    Raises:
        ValueError: ``dt_ps`` や単位・原点モードが不正、選択が空、フレームが2未満か等間隔でない、速度がないか非有限、lag 上限が範囲外、Cv(0) が0、選択式が不正か、一致する原子がない。

    Note:
        速度は保存されたまま使い、平均も重心ドリフトも引かない。拡散係数・粘性・スペクトルの推定はしない。
    """
    if not np.isfinite(dt_ps) or dt_ps <= 0:
        raise ValueError('dt_ps must be positive ps per LAMMPS timestep')
    if velocity_unit not in ('A/fs', 'A/ps'):
        raise ValueError("velocity_unit must be exactly 'A/fs' or 'A/ps'")
    if origin_mode not in ('multiple', 'single'):
        raise ValueError("origin_mode must be 'multiple' (Colab) or 'single' (frame 0 only)")
    if not selections:
        raise ValueError('Specify at least one atom selection')
    frames = _frames(trajectory, start_step=start_step, stop_step=stop_step, stride=stride)
    steps = _steps(trajectory)[frames]
    if len(frames) < 2 or not np.all(np.diff(steps) == steps[1] - steps[0]):
        raise ValueError('VACF needs at least two uniformly spaced selected frames')
    velocities = _velocities(trajectory, frames)
    if velocities is None:
        raise ValueError('VACF requires dump velocities (vx/vy/vz)')
    if not np.isfinite(velocities).all():
        raise ValueError('Non-finite velocities')
    if velocity_unit == 'A/fs':
        velocities = velocities * 1000.0
    frame_interval_ps = float(steps[1] - steps[0]) * dt_ps
    lags = np.arange(len(frames))
    lag_ps = lags * frame_interval_ps
    if max_lag_ps is not None:
        if not np.isfinite(max_lag_ps) or max_lag_ps <= 0:
            raise ValueError('max_lag_ps must be positive and finite')
        tolerance = 1e-9 * max(1.0, float(lag_ps[-1]))
        if max_lag_ps > lag_ps[-1] + tolerance:
            raise ValueError(f'max_lag_ps exceeds the available lag range 0..{lag_ps[-1]:g} ps')
        keep = lag_ps <= max_lag_ps + tolerance
        lags, lag_ps = lags[keep], lag_ps[keep]
    mean_velocity = velocities.mean(axis=1)
    mean_velocity_rms = float(np.sqrt(np.mean(np.sum(mean_velocity**2, axis=1))))
    table = {'lag_ps': lag_ps}
    info = []
    for index, expression in enumerate(selections):
        mask = _select(trajectory, expression)
        correlations, counts = _correlations(velocities[:, mask], lags, origin_mode)
        total = correlations.sum(axis=1)
        if not total[0] > 0:
            raise ValueError('Cv(0) is zero; normalized VACF is undefined')
        key = f'VACF_{index}'
        table[f'{key}_x_A2_ps2'] = correlations[:, 0]
        table[f'{key}_y_A2_ps2'] = correlations[:, 1]
        table[f'{key}_z_A2_ps2'] = correlations[:, 2]
        table[f'{key}_A2_ps2'] = total
        table[f'{key}_normalized'] = total / total[0]
        table[f'origin_atom_pairs_{index}'] = counts
        info.append({'selection': expression, 'column': f'{key}_A2_ps2',
                     'normalized': f'{key}_normalized', 'atoms': int(mask.sum())})
    limitations = ['Overlapping time origins are correlated; samples per lag decrease with lag',
                   'No equilibrium claim; no diffusion coefficient, viscosity or spectrum inferred',
                   'No mean or COM-drift subtraction; velocities are used as saved',
                   'mean_velocity_rms_A_ps is the RMS of the arithmetic mean velocity over all atoms '
                   '(equal-mass COM proxy only)']
    if origin_mode == 'single':
        limitations.insert(0, 'Single origin (frame 0) only; not the Colab multiple-origin average')
    return Result({'vacf': _frame(table)}, {
        'kind': 'vacf', 'title': 'Velocity autocorrelation function', 'method': 'custom',
        'origin_mode': origin_mode, 'velocity_unit_input': velocity_unit,
        'velocity_unit_output': 'A/ps', 'dt_ps': dt_ps, 'frame_interval_ps': frame_interval_ps,
        'source_dump': _source_dump(trajectory),
        'selected_steps': [int(steps[0]), int(steps[-1])], 'selected_frames': len(frames),
        'frame_stride': stride, 'mean_velocity_rms_A_ps': mean_velocity_rms,
        'drift_removed': False, 'vacf_selections': info, 'limitations': limitations,
    })
