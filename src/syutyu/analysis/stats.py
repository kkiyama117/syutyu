"""thermo の記述統計。平衡に達したとは自動で主張しない。
"""
import numpy as np
import polars as pl


def summarize_thermo(table, *, start_ps=None, stop_ps=None) -> pl.DataFrame:
    """thermo の区間について、平均・標準偏差・点数を返す。

    Args:
        table (polars.DataFrame): ``read_thermo`` が返す表。``time_ps`` 列が必要。
        start_ps (float | None): 集計する時間範囲（両端を含む）。省略時は全体。
        stop_ps (float | None): 集計する時間範囲（両端を含む）。省略時は全体。

    Returns:
        polars.DataFrame: ``time_ps`` と ``Step`` 以外の列ごとに1行。列は ``column``、``mean``、``std``（母集団標準偏差）、``samples``。

    Raises:
        TypeError: ``table`` が polars の DataFrame でない。
        ValueError: ``time_ps`` が非有限か減少、区間が不正、該当サンプルがない、列の長さや値が不正。

    Note:
        ``std`` は母集団標準偏差であり、平均の標準誤差ではない。記録されたサンプルは相関する。
    """
    if not isinstance(table, pl.DataFrame):
        raise TypeError('table must be a polars DataFrame')
    time = table['time_ps'].to_numpy().astype(float)
    if time.ndim != 1 or not len(time) or not np.isfinite(time).all() or np.any(np.diff(time) <= 0):
        raise ValueError('Thermo time_ps must be finite and strictly increasing')
    lower = time[0] if start_ps is None else start_ps
    upper = time[-1] if stop_ps is None else stop_ps
    if not np.isfinite([lower, upper]).all() or lower > upper:
        raise ValueError('Invalid thermo interval')
    selected = (time >= lower) & (time <= upper)
    if not selected.any():
        raise ValueError('No thermo samples in selected interval')
    records = []
    for name in table.columns:
        if name in ('time_ps', 'Step'):
            continue
        values = table[name].to_numpy().astype(float)
        if values.shape != time.shape or not np.isfinite(values).all():
            raise ValueError(f'Invalid thermo column: {name}')
        values = values[selected]
        records.append({'column': name, 'mean': float(values.mean()),
                        'std': float(values.std()), 'samples': len(values)})
    return pl.DataFrame(records, schema={'column': pl.String, 'mean': pl.Float64,
                                         'std': pl.Float64, 'samples': pl.Int64})
