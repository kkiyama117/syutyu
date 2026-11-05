"""LAMMPS の thermo 出力を、ブロックごとに分けて読む。単位系は仮定しない。
"""
from pathlib import Path

import numpy as np


def read_thermo_blocks(path) -> list[dict[str, np.ndarray]]:
    """``Step`` で始まり ``Loop time`` で終わる表を、ブロックごとに分けて読む。

    Args:
        path (str | os.PathLike): LAMMPS の log または screen 出力。

    Returns:
        list[dict[str, numpy.ndarray]]: ブロックごとの辞書（列名→numpy 配列）のリスト。列名と step はブロック間で違ってよい。

    Raises:
        ValueError: 列名の重複、列数や数値の不一致、``Loop time`` のない未完了ブロック。

    Note:
        WARNING などの表でない行は読み飛ばす。値は保存された単位のままで、換算しない。
    """
    blocks = []
    headers = None
    rows = []
    with Path(path).open(encoding='utf-8', errors='replace') as stream:
        for line_number, line in enumerate(stream, 1):
            fields = line.split()
            if not fields:
                continue
            if fields[0] == 'Step':
                if headers is not None:
                    raise ValueError(f'Unfinished thermo block before line {line_number}')
                if len(set(fields)) != len(fields):
                    raise ValueError(f'Duplicate thermo column at line {line_number}')
                headers, rows = fields, []
            elif headers is not None and line.lstrip().startswith('Loop time'):
                values = np.asarray(rows, dtype=float).reshape((-1, len(headers)))
                blocks.append({name: values[:, index].copy() for index, name in enumerate(headers)})
                headers = None
            elif headers is not None:
                try:
                    float(fields[0])
                except ValueError:
                    continue
                if len(fields) != len(headers):
                    raise ValueError(f'Wrong thermo column count at line {line_number}')
                try:
                    rows.append([float(value) for value in fields])
                except ValueError as error:
                    raise ValueError(f'Invalid thermo number at line {line_number}') from error
    if headers is not None:
        raise ValueError('Unfinished thermo block: missing Loop time marker')
    return blocks


def read_thermo(path, *, block=-1, dt_ps) -> dict[str, np.ndarray]:
    """1つの thermo ブロックを選び、``time_ps`` を加えて返す。

    Args:
        path (str | os.PathLike): LAMMPS の log または screen 出力。
        block (int): 選ぶブロックの番号。負なら末尾から数える。省略時は最後のブロック。
        dt_ps (float): LAMMPS の 1 step 当たりの ps。``timestep`` コマンドの値そのものではない。

    Returns:
        dict[str, numpy.ndarray]: 列名→numpy 配列の辞書。``time_ps = Step * dt_ps`` を加える。

    Raises:
        ValueError: ``dt_ps`` が正でない、ブロックがない、選択ブロックが空・非有限・step が増加しない、``time_ps`` 列と衝突する。
        IndexError: ブロック番号が範囲外。

    Note:
        ブロックは結合しない。reset の前後を1つの時系列として扱わない。
    """
    if not np.isfinite(dt_ps) or dt_ps <= 0:
        raise ValueError('dt_ps must be finite and positive (ps per timestep)')
    if not isinstance(block, int) or isinstance(block, bool):
        raise ValueError('block must be an integer index')
    blocks = read_thermo_blocks(path)
    if not blocks:
        raise ValueError('No thermo blocks found')
    if not -len(blocks) <= block < len(blocks):
        raise IndexError(f'Thermo block {block} out of range for {len(blocks)} blocks')
    table = blocks[block]
    if not len(table['Step']):
        raise ValueError('Selected thermo block is empty')
    if any(not np.isfinite(values).all() for values in table.values()):
        raise ValueError('Selected thermo block contains non-finite data')
    if not np.all(np.diff(table['Step']) > 0):
        raise ValueError('Selected thermo steps must be strictly ascending')
    if 'time_ps' in table:
        raise ValueError('Thermo column time_ps conflicts with the derived time column')
    with np.errstate(over='ignore', invalid='ignore'):
        time_ps = table['Step'] * dt_ps
    if not np.isfinite(time_ps).all():
        raise ValueError('Derived time_ps contains non-finite data')
    return {'time_ps': time_ps, **table}
