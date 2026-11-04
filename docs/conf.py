"""Sphinx configuration for the GitHub Pages site.

Build from the repository root:

    uv run --locked --group docs sphinx-build -W docs docs/_build/html
"""
import os
from pathlib import Path
import re
import subprocess

project = 'syutyu'
author = 'kkiyama117'
language = 'ja'

extensions = ['myst_parser']
myst_heading_anchors = 3

html_theme = 'furo'
html_title = 'syutyu ドキュメント'
exclude_patterns = ['_build']

DOCS = Path(__file__).resolve().parent
REPO = DOCS.parent
BLOB = 'https://github.com/kkiyama117/syutyu/blob/main/'
TRACKED = set(subprocess.check_output(['git', 'ls-files'], cwd=REPO, text=True).splitlines())
TRACKED_DIRS = {str(parent) for name in TRACKED for parent in Path(name).parents}
LINK = re.compile(r'(!?\[[^\]\n]*\])\((<[^>\n]*>|[^)\s]*)\)')


def outside_link(match, page):
    """Point links that leave docs/ at GitHub; drop links to untracked files.

    Sphinx would otherwise copy the target into the site, which would publish
    local records.
    """
    label, raw = match.groups()
    target = raw[1:-1] if raw.startswith('<') else raw
    if not target or re.match(r'^[a-zA-Z][\w+.-]*:|^[#/]', target):
        return match.group(0)
    location, _, anchor = target.partition('#')
    path = Path(os.path.normpath(page.parent / location))
    if path.is_relative_to(DOCS):
        return match.group(0)
    relative = path.relative_to(REPO).as_posix()
    if relative not in TRACKED and relative not in TRACKED_DIRS:
        return label[1:-1] if label.startswith('[') else label[2:-1]
    return f'{label}(<{BLOB}{relative}{"#" + anchor if anchor else ""}>)'


def rewrite_outside_links(app, docname, source):
    page = DOCS / app.env.doc2path(docname, False)
    source[0] = LINK.sub(lambda match: outside_link(match, page), source[0])


def setup(app):
    app.connect('source-read', rewrite_outside_links)
