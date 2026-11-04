"""thermo の記述統計。平衡に達したとは自動で主張しない。
"""
import numpy as np


def summarize_thermo(table, *, start_ps=None, stop_ps=None):
    """thermo の区間について、平均・標準偏差・点数を返す。

    引数:
        - `table` — `read_thermo` が返す辞書。`time_ps` 列が必要。
        - `start_ps`、`stop_ps` — 集計する時間範囲（両端を含む）。省略時は全体。
    戻り値:
        列名→`{'mean', 'std', 'samples'}` の辞書。`time_ps` と `Step` は含めない。
    例外:
        `ValueError` — `time_ps` が非有限か減少、区間が不正、該当サンプルがない、列の長さや値が不正。
    注意:
        `std` は母集団標準偏差であり、平均の標準誤差ではない。記録されたサンプルは相関する。
    """
    time = np.asarray(table['time_ps'], dtype=float)
    if time.ndim != 1 or not len(time) or not np.isfinite(time).all() or np.any(np.diff(time) <= 0):
        raise ValueError('Thermo time_ps must be finite and strictly increasing')
    lower = time[0] if start_ps is None else start_ps
    upper = time[-1] if stop_ps is None else stop_ps
    if not np.isfinite([lower, upper]).all() or lower > upper:
        raise ValueError('Invalid thermo interval')
    selected = (time >= lower) & (time <= upper)
    if not selected.any():
        raise ValueError('No thermo samples in selected interval')
    stats = {}
    for key, values in table.items():
        if key in ('time_ps', 'Step'):
            continue
        values = np.asarray(values, dtype=float)
        if values.shape != time.shape or not np.isfinite(values).all():
            raise ValueError(f'Invalid thermo column: {key}')
        values = values[selected]
        stats[key] = {'mean': float(values.mean()), 'std': float(values.std()), 'samples': len(values)}
    return stats
