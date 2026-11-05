"""Result の表を共通に描く。pyplot を import しない headless 描画。
"""
from pathlib import Path
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg

VIEWS = {'rdf': 'RDF', 'cn': 'Coordination', 'msd': 'MSD', 'vacf': 'VACF',
         'vacf-normalized': 'Normalized VACF', 'temperature': 'Temperature'}


def available_views(result):
    """Result から描ける表示の一覧を返す。

    Args:
        result (Result): ``Result``。

    Returns:
        dict[str, str]: キー→表示名の辞書。``rdf``、``cn``、``msd``、``vacf``、``vacf-normalized``、``temperature`` のうち描けるもの。

    Raises:
        ValueError: 描ける表がない。

    Note:
        ``temperature`` は ``thermo`` 表に ``Temp`` 列があるときだけ加わる。
    """
    keys = []
    if 'rdf' in result.tables:
        keys += ['rdf', 'cn']
    if 'msd' in result.tables:
        keys += ['msd']
    if 'vacf' in result.tables:
        keys += ['vacf', 'vacf-normalized']
    if 'thermo' in result.tables and 'Temp' in result.tables['thermo']:
        keys += ['temperature']
    if not keys:
        raise ValueError('No plottable RDF, MSD or temperature table')
    return {key: VIEWS[key] for key in keys}


def decorate_figure(fig, result):
    """図の題名・条件・注意書きを summary から書く。

    Args:
        fig (matplotlib.figure.Figure): Matplotlib の ``Figure``。
        result (Result): ``Result``。``summary`` の ``title``、``atoms``、``selected_steps``、``method``、``coordinate_frame`` などを読む。

    Note:
        可変セルの実験室座標 MSD では、補正済みの拡散係数ではないと明示する。
    """
    s = result.summary
    fig.patch.set_facecolor('#f5f7fa')
    fig.text(.05, .94, s.get('title', 'MD analysis'), fontsize=21, weight='bold', color='#18354a')
    details = []
    if 'atoms' in s:
        details.append(f'N = {s["atoms"]:,}')
    if 'selected_steps' in s:
        details.append(f'Steps {s["selected_steps"][0]}–{s["selected_steps"][1]} | {s["selected_frames"]} frames')
    if 'method' in s:
        details.append(f'method = {s["method"]}')
    fig.text(.05, .885, '  |  '.join(details), fontsize=11)
    warning = 'Equilibrium not confirmed; interpret using the actual simulation conditions.'
    if s.get('coordinate_frame', '').startswith('laboratory'):
        warning = 'Laboratory-frame MSD includes box motion — not a corrected diffusion estimate.'
    fig.text(.05, .84, warning, color='#9a3d15', fontsize=10)


def draw_view(ax, result, view=None, compare=True):
    """1つの軸に表示を描き、図の下へ出す注記を返す。

    Args:
        ax (matplotlib.axes.Axes): 描画先の Matplotlib の軸。内容は消してから描く。
        result (Result): ``Result``。
        view (str | None): 表示のキー。省略時は最初の表示。
        compare (bool): 参照線や過去フレームとの比較を重ねる。

    Returns:
        str: 図に添える注記の文字列。規格化や原点の前提を短く書く。

    Raises:
        ValueError: ``view`` が ``available_views`` にない。
    """
    views = available_views(result)
    view = view or next(iter(views))
    if view not in views:
        raise ValueError(f'View {view} is unavailable; choose {list(views)}')
    ax.clear()
    s = result.summary
    if view in ('rdf', 'cn'):
        data = result.tables['rdf']
        for pair in s['rdf_pairs']:
            label = f'{pair["center"]} → {pair["neighbor"]}'
            x, y = (data['r_center_A'], data[pair['g']]) if view == 'rdf' else (data['r_upper_A'], data[pair['cn']])
            ax.plot(x, y, lw=2, label=label)
        if view == 'rdf':
            if compare:
                ax.axhline(1, color='grey', ls=':', label='g = 1 reference')
            ax.set(xlabel='Distance r (Å)', ylabel='g(r)', title='Radial distribution function')
            hint = (f'Bin width = {s["bin_width_A"]:g} Å. Compare: g=1 and earlier frames if present.\n'
                    f'Normalization: {s.get("normalization", "recorded in the LAMMPS input")}.')
        else:
            if compare:
                for pair in s['rdf_pairs']:
                    ax.plot(data['r_upper_A'][-1], data[pair['cn']][-1], 'o')
            ax.set(xlabel='Upper bin edge r (Å)', ylabel='Neighbors per center atom', title='Cumulative coordination number')
            hint = 'Counted directly per center atom; chosen radii are not first-shell boundaries.\nFor changing volume, do not multiply average g(r) by one average density.'
    elif view == 'msd':
        data = result.tables['msd']
        for group in s['msd_selections']:
            ax.plot(data['lag_ps'], data[group['column']], lw=2, label=group['selection'])
        if compare:
            ax.axhline(0, color='grey', ls=':', label='MSD = 0')
        ax.set(xlabel='Lag time (ps)', ylabel='MSD (Å²)', title='Atomic mean squared displacement')
        if s.get('com_drift_removed'):
            note = ('Atomic MSD with system COM drift removed, not molecular COM; '
                    'inspect origin/averaging settings before diffusion fits.')
        else:
            note = ('Atomic MSD, not molecular COM; '
                    'inspect origin/averaging settings before diffusion fits.')
        hint = f'Unwrapping: {s.get("unwrap_source", "minimum-image increments (conditional)")}.\n{note}'
    elif view in ('vacf', 'vacf-normalized'):
        data = result.tables['vacf']
        normalized = view == 'vacf-normalized'
        for group in s['vacf_selections']:
            column = group['normalized'] if normalized else group['column']
            ax.plot(data['lag_ps'], data[column], lw=2, label=group['selection'])
        if compare:
            ax.axhline(0, color='grey', ls=':')
        ax.set(xlabel='Lag time (ps)', ylabel='Cᵥ(t) / Cᵥ(0)' if normalized else 'Cᵥ(t) (Å²/ps²)',
               title='Normalized velocity autocorrelation' if normalized else 'Velocity autocorrelation')
        hint = (f'Origins: {s.get("origin_mode", "specified in metadata")}; velocities used as saved.\n'
                'Single origin (LAMMPS) and multiple origins (Colab) are different estimators; no drift removal.')
    else:
        data = result.tables['thermo']
        ax.plot(data['time_ps'], data['Temp'], 'o-', label='Recorded temperature')
        if compare:
            if 'target_temperature_K' in s:
                ax.axhline(s['target_temperature_K'], ls='--', label=f'Target {s["target_temperature_K"]:g} K')
            else:
                ax.axhline(data['Temp'].mean(), ls=':', label='Sample mean (not target)')
        ax.set(xlabel='Time within selected thermo block (ps)', ylabel='Temperature (K)', title='Temperature history')
        hint = 'Time from explicit timestep; blocks separated at each thermo header.\nTemperature control alone does not establish equilibrium.'
    ax.grid(alpha=.18)
    ax.legend(fontsize=9)
    ax.margins(x=.02)
    return hint


def make_figure(result, view=None, compare=True):
    """描画済みの Figure を作って返す。保存はしない。

    Args:
        result (Result): ``Result``。
        view (str | None): 表示のキー。省略時は最初の表示。
        compare (bool): 参照線や過去フレームとの比較を重ねる。

    Returns:
        matplotlib.figure.Figure: Matplotlib の ``Figure``。10×6.5 インチ。

    Raises:
        ValueError: 描ける表がない、または ``view`` が使えない。
    """
    fig = Figure(figsize=(10, 6.5))
    FigureCanvasAgg(fig)
    ax = fig.subplots()
    fig.subplots_adjust(left=.12, right=.96, bottom=.25, top=.77)
    decorate_figure(fig, result)
    fig.text(.12, .10, draw_view(ax, result, view, compare), fontsize=10, color='#435567')
    return fig


def save_figure(result, output, view=None, compare=True, *, overwrite=False):
    """Result を PNG に描いて保存する。

    Args:
        result (Result): ``Result``。
        output (str | os.PathLike): 保存先のパス。親フォルダは作る。
        view (str | None): ``available_views`` のキー。省略時は最初の表示。
        compare (bool): 参照線や過去フレームとの比較を重ねる。
        overwrite (bool): 既存ファイルへの上書きを明示する。

    Returns:
        pathlib.Path: 保存先の ``Path``。解像度は 160 dpi。

    Raises:
        FileExistsError: 既存ファイルがあり ``overwrite`` を指定していない。
        ValueError: 描ける表がない、または ``view`` が使えない。

    Note:
        GUI を開かない。
    """
    output = Path(output)
    if output.exists() and not overwrite:
        raise FileExistsError(f'Refusing to overwrite {output}')
    fig = make_figure(result, view, compare)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=160)
    return output
