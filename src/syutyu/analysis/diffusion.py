"""明示した時間窓での3次元アインシュタイン拡散の見積もり。
"""
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class DiffusionFit:
    """1つの時間窓で MSD を直線近似した結果。

    属性:
        - `window_ps` — 近似に使った経過時間の範囲（両端を含む）。ps。
        - `samples` — 窓に入った点の数。
        - `slope_A2_ps` — MSD の傾き。Å²/ps。
        - `intercept_A2` — 切片。Å²。
        - `D_m2_s` — `slope / 6 × 1e-8` で求めた拡散係数。m²/s。
        - `r_squared` — 窓内での決定係数。線形性の目安であり、平衡や収束の証拠ではない。
        - `limitation` — この見積もりで主張しないことの既定の文。
    """
    window_ps: tuple[float, float]
    samples: int
    slope_A2_ps: float
    intercept_A2: float
    D_m2_s: float
    r_squared: float
    limitation: str = ('A linear fit is not evidence of equilibrium or convergence; '
                       'correlated MSD samples do not give independent error bars.')

    def predict(self, time_ps):
        """近似直線の値を返す。

        引数:
            - `time_ps` — 経過時間。ps。スカラーか配列。
        戻り値:
            `slope_A2_ps * time_ps + intercept_A2`。Å²。
        """
        return self.slope_A2_ps * np.asarray(time_ps) + self.intercept_A2


def fit_diffusion(time_ps, msd_A2, *, window_ps):
    """指定した時間窓で MSD を直線近似し、拡散係数を返す。

    引数:
        - `time_ps` — MSD 原点からの経過時間。ps。狭義単調増加、非負。
        - `msd_A2` — 同じ長さの MSD。Å²。非負。
        - `window_ps` — 近似に使う範囲 `(下限, 上限)`。ps。データの範囲内で順序どおりであること。
    戻り値:
        `DiffusionFit`。`D_m2_s` は `slope / 6 × 1e-8`。
    例外:
        `ValueError` — 長さ・形・有限性・単調性が条件を満たさない、窓が範囲外か逆順、窓内の点が3未満、傾きが負。
    注意:
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
    slope, intercept = np.polyfit(time[selected], msd[selected], 1)
    if slope < 0:
        raise ValueError('Negative fitted slope: this window cannot provide a physical diffusion estimate')
    residual = msd[selected] - (slope * time[selected] + intercept)
    total = np.sum((msd[selected] - msd[selected].mean()) ** 2)
    r_squared = float(1 - np.sum(residual ** 2) / total) if total else 1.0
    return DiffusionFit((float(lower), float(upper)), int(selected.sum()),
                        float(slope), float(intercept), float(slope / 6 * 1e-8), r_squared)
