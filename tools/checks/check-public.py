"""Fail when a tracked file of this public repository mentions the private course.

Run from the library repository root:
    uv run --locked tools/checks/check-public.py

The course materials, assignments, and course recipes live in a private
repository. This check keeps their paths, names, and the lecture host out of
every tracked text file here.
"""
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
SELF = Path(__file__).resolve().relative_to(ROOT).as_posix()
PATTERNS = [
    r'course/', r'assignments/', r'syutyu_course', r'argon[-_]course', r'argon[-_]slide',
    r'kadai', r'\.secret', r'amorphous\.tf', r'chiba-u', r'Classroom',
    r'講義', r'先進エネルギー', r'集中講義', r'授業', r'課題',
]
FOUND = re.compile('|'.join(PATTERNS), flags=re.IGNORECASE)


def tracked():
    listed = subprocess.check_output(['git', 'ls-files', '-z', '--cached', '--others', '--exclude-standard'],
                                     cwd=ROOT).split(b'\0')
    return [name.decode() for name in listed if name]


def main():
    hits = []
    for name in tracked():
        path = ROOT / name
        if name == SELF or not path.is_file():
            continue
        try:
            text = path.read_text(encoding='utf-8')
        except UnicodeDecodeError:
            continue
        hits += [f'{name}:{number}: {match.group(0)}'
                 for number, line in enumerate(text.splitlines(), 1)
                 for match in [FOUND.search(line)] if match]
    if hits:
        print('\n'.join(hits), file=sys.stderr)
        print(f'PRIVATE_REFERENCE {len(hits)} line(s)', file=sys.stderr)
        return 1
    print('public check ok')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
