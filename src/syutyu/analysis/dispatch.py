"""LAMMPS と自前計算を明示的に選ぶ。別の推定量へ黙って切り替えない。
"""
from collections.abc import Mapping

from ..io.dump import Trajectory, read_trajectory
from ..io.native import read_msd, read_rdf, read_vacf
from .msd import mean_squared_displacement
from .rdf import radial_distribution
from .vacf import velocity_autocorrelation


def _method(method):
    if method not in ('lammps', 'custom'):
        raise ValueError('method must be explicitly lammps or custom')


def _trajectory(source):
    return source if isinstance(source, Trajectory) else read_trajectory(source)


def rdf(source, *, method, normalization='finite', **options):
    """軌跡または LAMMPS の RDF 表から g(r) と積算配位数を返す。

    引数:
        - `source` — `method='custom'` なら `Trajectory` か dump のパス。`method='lammps'` なら `.gr` のパス。
        - `method` — `custom`（自前計算）か `lammps`（compute rdf の出力を読む）。必須。
        - `normalization` — `finite` か `population`。`lammps` では `finite` のみ。
        - `options` — `custom` では `pairs`・`rmax`・`bins`・`start_step`・`stop_step`・`stride`・`intermolecular`。`lammps` では `pairs`・`start_step`・`stop_step`。
    戻り値:
        `Result`。`tables['rdf']` に `r_center_A`、`r_upper_A`、`g_<i>`、`CN_<i>`。
    例外:
        `ValueError` — `method` が2つ以外、`lammps` で `population` を指定、その方法で使えないオプションを渡した。
        `TypeError` — `custom` で `pairs` や `rmax` を欠いた。
    注意:
        この入口の既定は有限N規格化（Colab と LAMMPS の流儀）。下位の `radial_distribution` は `normalization` を必須にしており、既定を継承しない。
    """
    _method(method)
    if method == 'custom':
        return radial_distribution(_trajectory(source), normalization=normalization, **options)
    if normalization != 'finite':
        raise ValueError('LAMMPS compute rdf uses finite normalization; population is custom-only')
    invalid = set(options) - {'pairs', 'start_step', 'stop_step'}
    if invalid:
        raise ValueError(f'Custom-only RDF options for lammps method: {sorted(invalid)}')
    result = read_rdf(source, **options)
    result.summary.update(method='lammps', rdf_normalization='finite',
                          normalization='LAMMPS compute rdf: N_A*N_B-N_overlap')
    return result


def msd(source, *, method, dt_ps, origin_mode=None, origin_step=None, **options):
    """自前計算の原子 MSD か、LAMMPS compute msd の出力を選んで返す。

    引数:
        - `source` — `custom` なら `Trajectory` か dump のパス。`lammps` ならファイル名か系列名→ファイル名の辞書。
        - `method` — `custom` か `lammps`。必須。
        - `dt_ps` — LAMMPS の 1 step 当たりの ps。
        - `origin_mode` — `multiple`（全開始時刻を平均）か `single`（最初のフレームを原点）。`lammps` は単一原点のみ。
        - `origin_step` — `lammps` では必須。`custom` では指定できない。
        - `options` — `custom` では `selections`（既定 `['all']`）、`start_step`、`stop_step`、`stride`、`laboratory`。`lammps` では `column`。
    戻り値:
        `Result`。`tables['msd']` に `lag_ps`、`MSD_<i>_A2`、`origin_atom_pairs_<i>`。
    例外:
        `ValueError` — `method` が2つ以外、`lammps` で `origin_step` がないか `multiple` を要求、`custom` で `origin_step` を渡した、使えないオプションを渡した。
    注意:
        原子の MSD であり、分子重心の MSD ではない。重なった時間原点は独立でなく、拡散係数は推定しない。
    """
    _method(method)
    if method == 'custom':
        if origin_step is not None:
            raise ValueError('For custom MSD select the reference frame using start_step, not origin_step')
        selections = options.pop('selections', ['all'])
        return mean_squared_displacement(_trajectory(source), selections, dt_ps=dt_ps,
                                         origin_mode=origin_mode or 'multiple', **options)
    if origin_mode not in (None, 'single'):
        raise ValueError('LAMMPS compute msd has a single origin; multiple is custom-only')
    if origin_step is None:
        raise ValueError('LAMMPS MSD requires explicit origin_step')
    invalid = set(options) - {'column'}
    if invalid:
        raise ValueError(f'Custom-only MSD options for lammps method: {sorted(invalid)}')
    files = source if isinstance(source, Mapping) else {'LAMMPS compute group': source}
    return read_msd(files, dt_ps=dt_ps, origin_step=origin_step, **options)


def vacf(source, *, method, dt_ps, velocity_unit, origin_mode=None, origin_step=None, **options):
    """自前計算の VACF か、LAMMPS compute vacf の出力を選んで返す。

    引数:
        - `source` — `custom` なら速度を含む `Trajectory` か dump のパス。`lammps` なら fix ave/time の表。
        - `method` — `custom` か `lammps`。必須。
        - `dt_ps` — LAMMPS の 1 step 当たりの ps。
        - `velocity_unit` — `A/fs` か `A/ps`。`A/fs` は二乗する前に `A/ps` へ換算する。
        - `origin_mode` — `multiple`（Colab の定義）か `single`（フレーム0のみ）。`lammps` は単一原点のみ。
        - `origin_step` — `lammps` では必須。`custom` では指定できない。
        - `options` — `custom` では `selections`（既定 `['all']`）、`start_step`、`stride`、`max_lag_ps`。`lammps` では `max_lag_ps`。
    戻り値:
        `Result`。`tables['vacf']` に `lag_ps`、成分 `VACF_<i>_x_A2_ps2` など、`VACF_<i>_A2_ps2`、`VACF_<i>_normalized`、`origin_atom_pairs_<i>`。
    例外:
        `ValueError` — 単位・原点モード・方法の指定が不正、`lammps` で `origin_step` がない、使えないオプションを渡した。
    注意:
        単一原点と複数原点は別の推定量であり、同じ値として比較しない。速度の平均や重心ドリフトの除去はしない。
    """
    _method(method)
    if method == 'custom':
        if origin_step is not None:
            raise ValueError('For custom VACF select the reference frame using start_step, not origin_step')
        selections = options.pop('selections', ['all'])
        return velocity_autocorrelation(_trajectory(source), selections, dt_ps=dt_ps,
                                        velocity_unit=velocity_unit, origin_mode=origin_mode or 'multiple', **options)
    if origin_mode not in (None, 'single'):
        raise ValueError('LAMMPS compute vacf has a single origin; multiple is custom-only')
    if origin_step is None:
        raise ValueError('LAMMPS VACF requires explicit origin_step')
    invalid = set(options) - {'max_lag_ps'}
    if invalid:
        raise ValueError(f'Custom-only VACF options for lammps method: {sorted(invalid)}')
    return read_vacf(source, dt_ps=dt_ps, origin_step=origin_step, velocity_unit=velocity_unit, **options)
