"""LAMMPS の compute 出力（rdf・msd・vacf）を、再計算せずに読む。
"""
from pathlib import Path

import numpy as np

from .store import Result


def read_rdf(path, *, pairs=None, start_step=None, stop_step=None):
    """.gr 形式の RDF を読み、選んだブロックの平均を等重みで返す。

    Args:
        path (str | os.PathLike): ``.gr`` ファイル。step 行数 の見出しの後に ``index r g_0 CN_0 [g_1 CN_1 ...]`` が続く。
        pairs (Sequence[tuple[str, str]] | None): 表示用の (中心, 隣接) ラベル。省略時は ``center_<i>`` と ``neighbor_<i>``。
        start_step (int | None): 平均に使うブロックの step 範囲（両端を含む）。
        stop_step (int | None): 平均に使うブロックの step 範囲（両端を含む）。

    Returns:
        Result: ``tables['rdf']`` は ``r_center_A``、``r_upper_A``、``g_<i>``、``CN_<i>``。

    Raises:
        ValueError: 見出し・行数・列数・step の単調増加・ビン中心の等間隔のいずれかが形式に合わない。

    Note:
        距離は Å とみなす。規格化と原子対の除外は LAMMPS の入力側の設定をそのまま残し、g を規格化し直したり、平均密度から CN を作り直したりしない。
        ``pairs`` は表示名であり、元素やタイプの対応を推測するものではない。
    """
    path = Path(path)
    records = []
    with path.open(encoding='utf-8') as stream:
        lines = ((number, line.split()) for number, line in enumerate(stream, 1)
                 if line.strip() and not line.lstrip().startswith('#'))
        for number, header in lines:
            try:
                if len(header) != 2:
                    raise ValueError
                step, count = (int(value) for value in header)
                if count < 1:
                    raise ValueError
            except ValueError as error:
                raise ValueError(f'Expected integer step/row-count header at line {number}') from error
            rows = []
            for index in range(count):
                try:
                    row_number, fields = next(lines)
                except StopIteration as error:
                    raise ValueError(f'Incomplete RDF block at step {step}') from error
                try:
                    row = [float(value) for value in fields]
                except ValueError as error:
                    raise ValueError(f'Invalid RDF number at line {row_number}') from error
                if (len(row) < 4 or len(row) % 2 or row[0] != index + 1
                        or (rows and len(row) != len(rows[0]))):
                    raise ValueError(f'Invalid RDF row at line {row_number}')
                rows.append(row)
            values = np.asarray(rows)
            if not np.isfinite(values).all() or np.any(values[:, 1:] < 0):
                raise ValueError(f'RDF values must be finite and nonnegative (step {step})')
            if records and step <= records[-1][0]:
                raise ValueError('RDF block steps must strictly increase')
            if records and (values.shape != records[0][1].shape
                            or not np.array_equal(values[:, 1], records[0][1][:, 1])):
                raise ValueError('RDF blocks must have identical bin centers and column counts')
            records.append((step, values))
    if not records:
        raise ValueError('No RDF blocks found')
    if start_step is not None and stop_step is not None and start_step > stop_step:
        raise ValueError('start_step must not exceed stop_step')
    selected = [(step, values) for step, values in records
                if (start_step is None or step >= start_step)
                and (stop_step is None or step <= stop_step)]
    if not selected:
        raise ValueError('No RDF blocks in selected interval')
    centers = selected[0][1][:, 1]
    width = 2 * centers[0]
    if width <= 0 or not np.allclose(centers, (np.arange(len(centers)) + .5) * width,
                                    rtol=1e-5, atol=1e-7):
        raise ValueError('Expected uniformly spaced RDF bin centers starting at half a bin width')
    average = np.mean([values[:, 2:] for _, values in selected], axis=0)
    count = average.shape[1] // 2
    labels = list(pairs) if pairs is not None else [(f'center_{i}', f'neighbor_{i}') for i in range(count)]
    if len(labels) != count or any(len(pair) != 2 for pair in labels):
        raise ValueError(f'Expected {count} center/neighbor label pairs')
    table = {'r_center_A': centers.copy(), 'r_upper_A': centers + width / 2}
    info = []
    for index, (center, neighbor) in enumerate(labels):
        g, cn = f'g_{index}', f'CN_{index}'
        table[g], table[cn] = average[:, 2 * index], average[:, 2 * index + 1]
        info.append({'center': center, 'neighbor': neighbor, 'g': g, 'cn': cn})
    return Result({'rdf': table}, {
        'kind': 'rdf', 'title': 'LAMMPS computed RDF', 'source_rdf': str(path.resolve()),
        'rdf_pairs': info, 'bin_width_A': float(width), 'selected_blocks': len(selected),
        'block_steps': [step for step, _ in selected],
        'averaging': 'Equal weights for selected fix ave/time block means',
        'limitations': ['Equilibrium not confirmed; blocks may be correlated',
                        'Normalization, exclusions and within-block averaging follow the LAMMPS input',
                        'Assumes distance unit Å; block means are not raw independent samples'],
    })


def read_msd(files, *, dt_ps, origin_step, column=-1):
    """compute msd の fix ave/time 表を、系列名ごとに読む。

    Args:
        files (Mapping[str, str | os.PathLike]): 系列名→ファイルの辞書。全系列で step の刻みが同じであること。
        dt_ps (float): LAMMPS の 1 step 当たりの ps。
        origin_step (int): compute の原点にした step。省略できない。
        column (int): MSD の列。step の列は選べない。

    Returns:
        Result: ``tables['msd']`` は ``lag_ps`` と ``MSD_<i>_A2``。

    Raises:
        ValueError: 列数や行数が足りない、非有限値、step が原点より前か増加しない、系列間で刻みが違う、step 列を選んだ。

    Note:
        選んだ列は Å² の全 MSD（成分や log ではない）を前提とする。座標の扱いと平均は LAMMPS の入力に従い、拡散係数や平衡性は推定しない。
    """
    if not files or not np.isfinite(dt_ps) or dt_ps <= 0:
        raise ValueError('Need named MSD files and a positive dt_ps')
    table, groups, reference_steps = {}, [], None
    for index, (name, path) in enumerate(files.items()):
        values = np.loadtxt(path, comments='#', ndmin=2)
        if values.shape[1] < 2 or not np.isfinite(values).all() or values.shape[0] < 1:
            raise ValueError('Expected finite numeric timestep + MSD rows')
        steps = values[:, 0]
        if np.any(np.diff(steps) <= 0) or np.any(steps < origin_step):
            raise ValueError('MSD steps must increase and follow the chosen origin')
        if not -values.shape[1] <= column < values.shape[1] or column % values.shape[1] == 0:
            raise ValueError('Select an MSD column, not the timestep column')
        if reference_steps is None:
            reference_steps = steps
            table['lag_ps'] = (steps - origin_step) * dt_ps
        if not np.array_equal(steps, reference_steps):
            raise ValueError('MSD files use different timestep grids; import separately')
        key = f'MSD_{index}_A2'
        table[key] = values[:, column]
        groups.append({'selection': name, 'column': key, 'source_name': Path(path).name})
    return Result({'msd': table}, {
        'kind': 'msd', 'title': 'LAMMPS computed MSD', 'msd_selections': groups,
        'method': 'lammps', 'origin_mode': 'single',
        'dt_ps': dt_ps, 'origin_step': origin_step, 'value_column': column,
        'unwrap_source': 'LAMMPS compute msd output',
        'limitations': ['Not a new multiple-origin estimate; values and averaging are defined by the input script',
                        'NPT cell motion and COM handling follow the original LAMMPS compute settings',
                        'No diffusion coefficient or equilibrium claim inferred'],
    })


def read_vacf(path, *, dt_ps, origin_step, velocity_unit, max_lag_ps=None):
    """compute vacf の表を読み、Å²/ps² へ換算して返す。

    Args:
        path (str | os.PathLike): step と vx*vx0、vy*vy0、vz*vz0、合計の5列を持つ表。
        dt_ps (float): LAMMPS の 1 step 当たりの ps。
        origin_step (int): compute の原点にした step。先頭の保存 step と一致すること。
        velocity_unit (str): 入力速度の単位。``A/fs`` か ``A/ps``。
        max_lag_ps (float | None): 残す lag の上限（両端を含む）。省略時は全 lag。

    Returns:
        Result: ``tables['vacf']`` は ``lag_ps``、成分 ``VACF_0_x_A2_ps2`` など、``VACF_0_A2_ps2``、``VACF_0_normalized``。

    Raises:
        ValueError: 単位が2つ以外、5列でない、step が原点から始まらないか増加しない、成分の和が合計と合わない、lag 上限が範囲外、Cv(0) が正でない。

    Note:
        単一原点の値であり、複数原点の平均ではない。速度の平均や重心ドリフトの除去はしない。
    """
    if not np.isfinite(dt_ps) or dt_ps <= 0:
        raise ValueError('dt_ps must be finite and positive')
    if velocity_unit not in ('A/fs', 'A/ps'):
        raise ValueError("velocity_unit must be A/fs or A/ps")
    if type(origin_step) is not int or origin_step < 0:
        raise ValueError('origin_step must be a nonnegative integer')
    values = np.loadtxt(path, comments='#', ndmin=2)
    if values.shape[1] != 5 or len(values) < 2 or not np.isfinite(values).all():
        raise ValueError('Expected finite timestep + four compute vacf columns, at least two rows')
    steps = values[:, 0]
    if (not np.array_equal(steps, np.rint(steps)) or np.any(np.diff(steps) <= 0)
            or steps[0] != origin_step):
        raise ValueError('VACF steps must increase and include origin_step as the first saved row')
    if not np.allclose(values[:, 1:4].sum(axis=1), values[:, 4], rtol=2e-5, atol=1e-12):
        raise ValueError('VACF components do not sum to the total; check column format')
    time = (steps - origin_step) * dt_ps
    keep = np.ones(len(time), dtype=bool)
    if max_lag_ps is not None:
        if not np.isfinite(max_lag_ps) or not 0 < max_lag_ps <= time[-1] + 1e-9:
            raise ValueError('max_lag_ps must lie within the available positive lag range')
        keep = time <= max_lag_ps + 1e-9
    correlations = values[keep, 1:] * (1e6 if velocity_unit == 'A/fs' else 1)
    if correlations[0, 3] <= 0:
        raise ValueError('Cv(0) must be positive for normalization')
    table = {'lag_ps': time[keep],
             **{f'VACF_0_{axis}_A2_ps2': correlations[:, index] for index, axis in enumerate('xyz')},
             'VACF_0_A2_ps2': correlations[:, 3],
             'VACF_0_normalized': correlations[:, 3] / correlations[0, 3]}
    return Result({'vacf': table}, {
        'kind': 'vacf', 'title': 'LAMMPS velocity autocorrelation', 'method': 'lammps',
        'origin_mode': 'single', 'origin_step': origin_step, 'dt_ps': dt_ps,
        'velocity_unit_input': velocity_unit, 'velocity_unit_output': 'A/ps',
        'source_vacf': str(Path(path).resolve()),
        'vacf_selections': [{'selection': 'LAMMPS compute group', 'column': 'VACF_0_A2_ps2',
                             'normalized': 'VACF_0_normalized'}],
        'limitations': ['Single compute origin only; not the Colab multiple-origin average',
                        'Coordinate units, compute group and within-output averaging follow the LAMMPS input',
                        'No equilibrium or convergence claim; no drift removal is performed by this reader'],
    })
