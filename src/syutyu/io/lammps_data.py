"""LAMMPS の atom_style charge 用 data ファイルを、分率座標から書く。
"""
from dataclasses import dataclass
import math
from pathlib import Path


@dataclass(frozen=True)
class Site:
    """直方晶の単位胞にある1原子。

    属性:
        - `type` — LAMMPS の原子タイプ。1 以上。
        - `charge` — 電荷。e。
        - `fractional` — 分率座標。各成分は [0, 1)。セル長を掛けると Å の座標になる。
    """
    type: int
    charge: float
    fractional: tuple[float, float, float]


def _float(value, message):
    try:
        return float(value)
    except (TypeError, ValueError):
        raise ValueError(message) from None


def _validated(lengths, sites, masses, title):
    """全入力を検査し、セル長・原子・タイプごとの質量を数値にして返す。"""
    if not isinstance(title, str) or not title.strip() or '\n' in title or '\r' in title:
        raise ValueError('title must be one non-empty line')
    try:
        cell = tuple(float(value) for value in lengths)
    except (TypeError, ValueError):
        raise ValueError('lengths must be three finite values > 0') from None
    if len(cell) != 3 or not all(math.isfinite(value) and value > 0 for value in cell):
        raise ValueError('lengths must be three finite values > 0')

    sites = list(sites)
    if not sites:
        raise ValueError('sites must not be empty')
    for index, site in enumerate(sites):
        if not isinstance(site.type, int) or site.type < 1:
            raise ValueError(f'site {index}: type must be an int >= 1')
        message = f'site {index}: fractional coordinates must be three values in [0, 1)'
        try:
            point = tuple(_float(value, message) for value in site.fractional)
        except TypeError:
            raise ValueError(message) from None
        if len(point) != 3 or not all(0.0 <= value < 1.0 for value in point):
            raise ValueError(message)
        if not math.isfinite(_float(site.charge, f'site {index}: charge must be finite')):
            raise ValueError(f'site {index}: charge must be finite')

    mass_by_type = {}
    for atom_type in sorted(masses):
        message = f'mass of type {atom_type} must be > 0'
        mass = _float(masses[atom_type], message)
        if not math.isfinite(mass) or mass <= 0:
            raise ValueError(message)
        mass_by_type[atom_type] = mass
    used = sorted({site.type for site in sites})
    missing = [atom_type for atom_type in used if atom_type not in mass_by_type]
    if missing:
        raise ValueError(f'types without a mass: {missing}')
    unused = [atom_type for atom_type in sorted(mass_by_type) if atom_type not in set(used)]
    if unused:
        raise ValueError(f'masses for unused types: {unused}')
    return cell, sites, mass_by_type


def _text(title, cell, sites, mass_by_type):
    """data ファイルの本文を、指定の書式と行順で組み立てる。"""
    types = sorted(mass_by_type)
    lines = [title, '', f'{len(sites)} atoms', f'{len(types)} atom types', '']
    lines += [f'0.0 {length:.10f} {axis}lo {axis}hi'
              for axis, length in zip('xyz', cell)]
    lines += ['', 'Masses', '']
    lines += [f'{atom_type} {mass_by_type[atom_type]:.6f}' for atom_type in types]
    lines += ['', 'Atoms # charge', '']
    for atom_id, site in enumerate(sites, 1):
        x, y, z = (float(value) * length
                   for value, length in zip(site.fractional, cell))
        lines.append(f'{atom_id} {site.type} {float(site.charge):+.6f}'
                     f' {x:.10f} {y:.10f} {z:.10f}')
    return '\n'.join(lines) + '\n'


def write_charge_data(path, lengths, sites, masses, *, title='syutyu charge data',
                      overwrite=False) -> Path:
    """分率座標の単位胞から LAMMPS の atom_style charge 用 data ファイルを書く。

    引数:
        - `path` — 書き出すファイル。親ディレクトリは存在すること。
        - `lengths` — 直方晶セルの (a, b, c)。Å。三つの有限な正の値。
        - `sites` — 空でない `Site` の列。原子 ID はこの順に 1 から振る。
        - `masses` — 原子タイプから質量（g/mol）への辞書。`sites` が使うタイプをすべて含み、使わないタイプを含まないこと。
        - `title` — 1行目の見出し。改行を含まない空でない1行。
        - `overwrite` — 既存ファイルを上書きすることを明示する。省略時は拒否する。
    戻り値:
        書き出したファイルの `Path`。
    例外:
        `ValueError` — 長さ・原子タイプ・分率座標・電荷・質量・タイプの対応・`title` のいずれかが規則に合わない。
        `FileExistsError` — `path` が存在し `overwrite` を指定していない。
        `FileNotFoundError` — `path` の親ディレクトリがない。
    注意:
        セルは原点 0 の直方晶で、角度は指定しない。電荷は e、長さは Å、質量は g/mol。
    """
    path = Path(path)
    cell, sites, mass_by_type = _validated(lengths, sites, masses, title)
    with path.open('w' if overwrite else 'x', encoding='utf-8', newline='\n') as stream:
        stream.write(_text(title, cell, sites, mass_by_type))
    return path
