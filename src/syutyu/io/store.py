"""全解析が共有する現行の結果保存形式。過去形式の読み替えは行わない。
"""
from dataclasses import dataclass
from pathlib import Path
import json

import polars as pl

SCHEMA = 3
TABLES = {'rdf', 'msd', 'vacf', 'thermo'}


def _frame(columns):
    """列名→配列の辞書を、系列の番号を持つ列を Int64、それ以外を Float64 にした表にする。"""
    schema = {name: (pl.Int64 if name.startswith(('ordered_counts_', 'origin_atom_pairs_'))
                     else pl.Float64) for name in columns}
    return pl.DataFrame(columns, schema=schema)


def _read_table(path):
    """表として読める CSV を返し、列も行もないファイルを拒否する。"""
    try:
        table = pl.read_csv(path)
    except pl.exceptions.NoDataError:
        table = pl.DataFrame()
    if not table.columns or table.height == 0:
        raise ValueError(f'Empty or unnamed CSV: {path}')
    return table


@dataclass
class Result:
    """解析結果の表と要約をまとめた値。

    Attributes:
        tables (dict[str, pl.DataFrame]): 表の名前（``rdf``、``msd``、``vacf``、``thermo``）から polars の DataFrame への辞書。
        summary (dict): 方法・単位・選択・限界を記録する辞書。保存時に ``schema_version`` と ``tables`` を加える。
    """
    tables: dict[str, pl.DataFrame]
    summary: dict


def save_result(result, output, *, overwrite=False):
    """Result を schema3 の CSV と summary.json として保存する。

    Args:
        result (Result): 保存する ``Result``。表は polars の DataFrame で、行があり、``rdf``・``msd``・``vacf``・``thermo`` のいずれかを含むこと。
        output (str | os.PathLike): 保存先ディレクトリ。なければ作る。
        overwrite (bool): 既存の非空ディレクトリへ書くことを明示する。省略時は拒否する。

    Returns:
        pathlib.Path: 保存先の ``Path``。

    Raises:
        ValueError: 表が空、表が polars の DataFrame でない、表に行がない、未知の表名を含む。summary に JSON へ書けない値（NaN など）がある。
        FileExistsError: ``output`` が非空で ``overwrite`` を指定していない。

    Note:
        書き込み前に既存ファイルを消さない。``overwrite=True`` でも同じ名前の CSV は上書きされる。
        ``summary.json`` は ``schema_version: 3`` と表の一覧を持ち、古い形式は読めない。
    """
    output = Path(output)
    if not result.tables or not set(result.tables) <= TABLES:
        raise ValueError('Result contains unknown or empty tables')
    if not all(isinstance(table, pl.DataFrame) for table in result.tables.values()):
        raise ValueError('Result tables must be polars DataFrames')
    if any(table.height == 0 for table in result.tables.values()):
        raise ValueError('Result tables must have rows')
    if output.exists() and any(output.iterdir()) and not overwrite:
        raise FileExistsError(f'Refusing to overwrite nonempty result directory: {output}')
    summary = {**result.summary, 'schema_version': SCHEMA, 'tables': list(result.tables)}
    serialized = json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + '\n'
    output.mkdir(parents=True, exist_ok=True)
    for name, table in result.tables.items():
        table.write_csv(output / f'{name}.csv')
    (output / 'summary.json').write_text(serialized)
    return output


def load_result(directory):
    """schema3 の保存ディレクトリを読み、Result を返す。

    Args:
        directory (str | os.PathLike): ``save_result`` が作ったディレクトリ。

    Returns:
        Result: ``tables`` の各表は polars の DataFrame。

    Raises:
        ValueError: ``schema_version`` が 3 でない、表の一覧が空、未知の表名を含む、CSV に列か行がない。
        FileNotFoundError: ``summary.json`` がない。

    Note:
        過去形式を読み替える処理はない。CSV は数値前提で、文字列の列は読めない。
    """
    directory = Path(directory)
    summary = json.loads((directory / 'summary.json').read_text())
    if summary.get('schema_version') != SCHEMA:
        raise ValueError('Unsupported result format; regenerate with the current CLI (no legacy adapter)')
    tables = summary.get('tables', [])
    if not tables or not set(tables) <= TABLES:
        raise ValueError('Unknown result tables')
    return Result({name: _read_table(directory / f'{name}.csv') for name in tables}, summary)
