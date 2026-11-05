"""明示した時間窓での3次元アインシュタイン拡散の見積もり。
"""
from dataclasses import dataclass

import numpy as np
from scipy import stats


@dataclass(frozen=True)
class DiffusionFit:
    """1つの時間窓で MSD を直線近似した結果。

    Attributes:
        window_ps (tuple[float, float]): 近似に使った経過時間の範囲（両端を含む）。ps。
        samples (int): 窓に入った点の数。
        regression (object): ``scipy.stats.linregress`` が返す ``LinregressResult``。
            ``slope`` は MSD の傾きで Å²/ps、``intercept`` は切片で Å²、``rvalue ** 2`` は決定係数。
            ``stderr``、``intercept_stderr``、``pvalue`` は標本が独立だと仮定した値であり、
            相関した MSD の点は独立ではない。
        D_m2_s (float): ``slope / 6 × 1e-8`` で求めた拡散係数。m²/s。
        limitation (str): この見積もりで主張しないことの既定の文。
    """
    window_ps: tuple[float, float]
    samples: int
    regression: object
    D_m2_s: float
    limitation: str = ('A linear fit is not evidence of equilibrium or convergence; '
                       'correlated MSD samples do not give independent error bars.')


def fit_diffusion(time_ps, msd_A2, *, window_ps):
    """指定した時間窓で MSD を直線近似し、拡散係数を返す。

    Args:
        time_ps (array_like): MSD 原点からの経過時間。ps。狭義単調増加、非負。
        msd_A2 (array_like): 同じ長さの MSD。Å²。非負。
        window_ps (tuple[float, float]): 近似に使う範囲 ``(下限, 上限)``。ps。データの範囲内で順序どおりであること。

    Returns:
        DiffusionFit: ``D_m2_s`` は ``slope / 6 × 1e-8``。

    Raises:
        ValueError: 長さ・形・有限性・単調性が条件を満たさない、窓が範囲外か逆順、窓内の点が3未満、傾きが負。

    Note:
        窓は物理的な拡散領域を自分で選ぶ。最小二乗は自動で窓を選ばず、平衡や収束、独立な誤差を主張しない。
    """
    time, msd = np.asarray(time_ps, dtype=float), np.asarray(msd_A2, dtype=float)
    if (time.ndim != 1 or msd.shape != time.shape or len(time) < 3
            or not np.isfinite(time).all() or not np.isfinite(msd).all()
            or np.any(time < 0) or np.any(msd < 0) or np.any(np.diff(time) <= 0)):
        raise ValueError('Need finite, nonnegative MSD and strictly increasing elapsed times')
    lower, upper = window_ps
    if not np.isfinite([lower, upper]).all() or not time[0] <= lower < upper <= time[-1]:
        raise ValueError('Fit window must be ordered and contained in the available elapsed times')
    selected = (time >= lower - 1e-10) & (time <= upper + 1e-10)
    if selected.sum() < 3:
        raise ValueError('Fit requires at least three samples within the window')
    regression = stats.linregress(time[selected], msd[selected])
    if regression.slope < 0:
        raise ValueError('Negative fitted slope: this window cannot provide a physical diffusion estimate')
    return DiffusionFit((float(lower), float(upper)), int(selected.sum()), regression,
                        float(regression.slope / 6 * 1e-8))
