"""Generate the API reference from the public docstrings.

Run from the repository root:

    uv run --locked tools/docs/api_reference.py --check
    uv run --locked tools/docs/api_reference.py --write

``--check`` fails when the committed reference no longer matches the code: a
renamed name, a changed signature, a new export, or an edited docstring.
``--write`` regenerates ``docs/reference/api.md``. The output holds no time, path
or machine information, so a clean tree regenerates the same bytes.

Docstring format parsed here (a program, not a model):

    One-line summary.

    引数:
        - ``name`` — type, unit, default, required or optional.
    戻り値:
        ...
    例外:
        ...
    注意:
        ...

A dataclass class docstring also carries ``属性:`` with one bullet per field;
the check fails when a field has no description. Keep a section label at the
left margin and indent its body.
"""
import argparse
import importlib
import inspect
import json
import re
import sys
import textwrap
from dataclasses import fields, is_dataclass
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
REFERENCE = REPO / 'docs/reference/api.md'
LABELS = ('引数', '戻り値', '例外', '注意', '属性')
SURFACE = (
    ('syutyu.io.store', ('Result', 'save_result', 'load_result')),
    ('syutyu.io.dump', ('Trajectory', 'read_trajectory', 'frame_indices')),
    ('syutyu.io.thermo', ('read_thermo', 'read_thermo_blocks')),
    ('syutyu.io.native', ('read_rdf', 'read_msd', 'read_vacf')),
    ('syutyu.io.lammps_data', ('Site', 'write_charge_data')),
    ('syutyu.analysis.dispatch', ('rdf', 'msd', 'vacf')),
    ('syutyu.analysis.diffusion', ('DiffusionFit', 'fit_diffusion')),
    ('syutyu.analysis.stats', ('summarize_thermo',)),
    ('syutyu.analysis.rdf', ('radial_distribution',)),
    ('syutyu.analysis.msd', ('mean_squared_displacement',)),
    ('syutyu.analysis.vacf', ('velocity_autocorrelation',)),
    ('syutyu.run', ('run_lammps',)),
    ('syutyu.report.figures', ('save_figure', 'available_views', 'draw_view',
                               'decorate_figure', 'make_figure')),
    ('syutyu.report.bundle', ('export_vmd', 'iter_dump_blocks')),
    ('syutyu.kernels', ('minimum_image', 'pair_rdf', 'unwrap', 'msd_multi_origin')),
)


def parse_docstring(obj):
    """Split one docstring into its summary paragraph and labelled sections."""
    text = inspect.getdoc(obj) or ''
    summary, sections, body = [], [], None
    for line in text.splitlines():
        stripped = line.strip()
        label = stripped[:-1].strip() if stripped.endswith(':') else ''
        if not line[:1].isspace() and label in LABELS:
            body = []
            sections.append((label, body))
        elif body is None:
            summary.append(line)
        else:
            body.append(line)
    cleaned = [(label, textwrap.dedent('\n'.join(lines)).strip())
               for label, lines in sections]
    return textwrap.dedent('\n'.join(summary)).strip(), cleaned


def public_members(cls):
    """Public methods and properties declared on the class itself."""
    members = []
    for name, value in vars(cls).items():
        if name.startswith('_'):
            continue
        if isinstance(value, property) and value.fget is not None:
            members.append((value.fget.__code__.co_firstlineno, name, value.fget, True))
        elif inspect.isfunction(value):
            members.append((value.__code__.co_firstlineno, name, value, False))
    return [(name, function, is_property)
            for _line, name, function, is_property in sorted(members)]


def check_callable(name, function, is_property):
    problems = []
    summary, sections = parse_docstring(function)
    if not summary:
        problems.append(f'{name}: 要約行がない')
    parameters = list(inspect.signature(function).parameters.values())
    if is_property:
        parameters = parameters[1:]
    if parameters and '引数' not in {label for label, _body in sections}:
        problems.append(f'{name}: 引数セクションがない')
    return problems


def check_entry(dotted, name, obj):
    problems = []
    summary, sections = parse_docstring(obj)
    if not summary:
        problems.append(f'{dotted}.{name}: 要約行がない')
    if inspect.isclass(obj):
        if not is_dataclass(obj):
            problems.append(f'{dotted}.{name}: 属性を列挙するため dataclass にする')
        else:
            documented = set(re.findall(r'^-\s+`([^`]+)`',
                                        dict(sections).get('属性', ''), flags=re.M))
            for field in fields(obj):
                if field.name not in documented:
                    problems.append(f'{dotted}.{name}: 属性 {field.name} の記述がない')
        for method_name, function, is_property in public_members(obj):
            problems += check_callable(f'{dotted}.{name}.{method_name}', function, is_property)
    else:
        problems += check_callable(f'{dotted}.{name}', obj, False)
    return problems


def collect():
    """Resolve every documented object and report docstring problems."""
    import syutyu
    facade = set(syutyu.__all__)
    entries, problems, seen = [], [], {}
    for dotted, names in SURFACE:
        module = importlib.import_module(dotted)
        for name in names:
            if not hasattr(module, name):
                problems.append(f'{dotted}: {name} が見つからない')
                continue
            if name in seen:
                problems.append(f'{name}: {seen[name]} と {dotted} で重複している')
            seen[name] = dotted
            entries.append((dotted, name, getattr(module, name), name in facade))
    for name in sorted(facade - set(seen)):
        problems.append(f'{name}: syutyu.__all__ にあるが SURFACE にない')
    for dotted, name, obj, _public in entries:
        problems += check_entry(dotted, name, obj)
    return entries, problems


def format_signature(function, name, drop_self=False):
    """Signature text without annotations; argument types live in 引数."""
    signature = inspect.signature(function)
    parameters = list(signature.parameters.values())
    if drop_self and parameters and parameters[0].name in ('self', 'cls'):
        parameters = parameters[1:]
    parameters = [item.replace(annotation=inspect.Parameter.empty)
                  for item in parameters]
    signature = signature.replace(parameters=parameters,
                                  return_annotation=inspect.Signature.empty)
    return f'{name}{signature}'


def anchor(name):
    return name.lower()


def render_section(label, body):
    """One labelled section; 例外 lists one entry per line."""
    if label == '例外':
        body = '\n'.join(f'- {line}' for line in body.splitlines() if line.strip())
    return [f'**{label}**', '', body, '']


def render_entry(dotted, name, obj, public):
    lines = [f'### {name}', '']
    if inspect.isclass(obj):
        lines += ['```python', f'class {name}', '```', '']
    else:
        lines += ['```python', format_signature(obj, name), '```', '']
    summary, sections = parse_docstring(obj)
    lines += [summary, '']
    if public:
        lines += [f'公開API: `from syutyu import {name}`', '']
    else:
        lines += [f'`from {dotted} import {name}`', '']
    for label, body in sections:
        lines += render_section(label, body)
    if inspect.isclass(obj):
        for method_name, function, is_property in public_members(obj):
            lines += [f'#### {method_name}', '']
            if is_property:
                lines += ['```python', method_name, '```', '']
            else:
                lines += ['```python',
                          format_signature(function, method_name, drop_self=True),
                          '```', '']
            member_summary, member_sections = parse_docstring(function)
            lines += [member_summary, '']
            for label, body in member_sections:
                lines += render_section(label, body)
    return lines


def render(entries):
    facade_count = sum(1 for _dotted, _name, _obj, public in entries if public)
    lines = [
        '# APIリファレンス',
        '',
        'このファイルは docstring から生成する。直接編集しない。',
        '',
        '```bash',
        'uv run --locked tools/docs/api_reference.py --check',
        '```',
        '',
        '- 対象は `syutyu` の公開クラスと関数。`cli.py` の内部、'
        '`syutyu.io.csvtable`、`_` で始まる名前は対象外。',
        f'- `from syutyu import ...` で使える{facade_count}名には「公開API」と記す。'
        'それ以外は submodule を明示して import する。',
        '- 単位・計算方法・実行例は [利用ガイド](../guide/usage.md)、'
        '[解析方法](../guide/analysis-methods.md)。',
        '',
        '## 一覧',
        '',
        '| 名前 | 種別 | 概要 |',
        '| --- | --- | --- |',
    ]
    for _dotted, name, obj, _public in entries:
        summary, _sections = parse_docstring(obj)
        kind = 'クラス' if inspect.isclass(obj) else '関数'
        lines.append(f'| [`{name}`](#{anchor(name)}) | {kind} | {summary} |')
    lines.append('')
    current = None
    for dotted, name, obj, public in entries:
        if dotted != current:
            current = dotted
            module = importlib.import_module(dotted)
            lines += [f'## {dotted}', '', (inspect.getdoc(module) or '').strip(), '']
        lines += render_entry(dotted, name, obj, public)
    return '\n'.join(lines).rstrip() + '\n'


def first_difference(expected, current):
    expected_lines, current_lines = expected.splitlines(), current.splitlines()
    for number, (left, right) in enumerate(zip(expected_lines, current_lines), 1):
        if left != right:
            return number
    if len(expected_lines) != len(current_lines):
        return min(len(expected_lines), len(current_lines)) + 1
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    action = parser.add_mutually_exclusive_group()
    action.add_argument('--write', action='store_true',
                        help='regenerate docs/reference/api.md')
    action.add_argument('--check', action='store_true',
                        help='verify the committed reference (default)')
    action.add_argument('--surface', action='store_true',
                        help='print the documented surface as JSON')
    parser.add_argument('--output', type=Path, help='optional verification receipt')
    args = parser.parse_args(argv)

    entries, problems = collect()
    if problems:
        for problem in problems:
            print(problem, file=sys.stderr)
        return 2
    if args.surface:
        print(json.dumps([{'module': dotted, 'name': name,
                           'kind': 'class' if inspect.isclass(obj) else 'function',
                           'public': public, 'signature': (
                               f'class {name}' if inspect.isclass(obj)
                               else format_signature(obj, name))}
                          for dotted, name, obj, public in entries],
                         ensure_ascii=False, indent=2))
        return 0

    text = render(entries)
    if args.write:
        REFERENCE.parent.mkdir(parents=True, exist_ok=True)
        REFERENCE.write_text(text, encoding='utf-8')
        status, difference = 'written', 0
    else:
        current = REFERENCE.read_text(encoding='utf-8') if REFERENCE.exists() else ''
        difference = first_difference(text, current)
        status = 'ok' if not difference else 'drift'

    receipt = {
        'entries': len(entries),
        'facade_entries': sum(1 for _d, _n, _o, public in entries if public),
        'dataclass_fields': sum(len(fields(obj)) for _d, _n, obj, _p in entries
                                if inspect.isclass(obj) and is_dataclass(obj)),
        'status': status,
        'first_difference_line': difference,
        'reference': str(REFERENCE.relative_to(REPO)),
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + '\n',
                               encoding='utf-8')
    print(json.dumps(receipt, ensure_ascii=False, indent=2))
    if status == 'drift':
        print(f'{receipt["reference"]} は現在の docstring と一致しない'
              f'（最初の違い: {difference} 行目）。--write で再生成する。', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
