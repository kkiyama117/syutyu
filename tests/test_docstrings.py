"""The public docstrings are Google style, and the reference lists every object."""
import importlib
import inspect
import pkgutil
import re
import unittest
from dataclasses import fields, is_dataclass
from pathlib import Path

import syutyu

SKIPPED = {'syutyu.cli', 'syutyu.io.csvtable'}
REFERENCE = Path(__file__).resolve().parents[1] / 'docs/reference/api.md'
CUSTOM_LABELS = ('引数:', '戻り値:', '例外:', '注意:', '属性:')
KANA = r'[\u3040-\u30fa\u30fc-\u30ff\u4e00-\u9fff]'
LONE_BACKTICK = re.compile(r'(?<!`)`(?!`)')
AFTER_KANA = re.compile(KANA + '``')
BEFORE_KANA = re.compile('``' + KANA)


def collect():
    """Return the public objects and the class members declared on them."""
    objects, members = {}, {}
    for found in pkgutil.walk_packages(syutyu.__path__, 'syutyu.'):
        if found.ispkg or found.name in SKIPPED:
            continue
        module = importlib.import_module(found.name)
        for name, value in vars(module).items():
            if name.startswith('_') or getattr(value, '__module__', None) != found.name:
                continue
            if not (inspect.isfunction(value) or inspect.isclass(value)):
                continue
            objects[f'{found.name}.{name}'] = (found.name, name, value)
            if inspect.isclass(value):
                for member, declared in vars(value).items():
                    if member.startswith('_'):
                        continue
                    if isinstance(declared, property) and declared.fget is not None:
                        members[f'{found.name}.{name}.{member}'] = declared.fget
                    elif inspect.isfunction(declared):
                        members[f'{found.name}.{name}.{member}'] = declared
    return objects, members


def documented(value):
    return inspect.getdoc(value) or ''


def section_lines(lines, header):
    start = next((index for index, line in enumerate(lines) if line.strip() == header), None)
    return [] if start is None else lines[start + 1:]


class PublicDocstringTest(unittest.TestCase):
    def test_public_objects_have_docstrings(self):
        objects, members = collect()
        self.assertGreaterEqual(len(objects) + len(members), 30)
        for name, value in {**{key: item[2] for key, item in objects.items()}, **members}.items():
            with self.subTest(name=name):
                self.assertTrue(documented(value).strip())

    def test_reference_lists_every_public_object(self):
        objects, _members = collect()
        expected = {
            f'syutyu.{name}' if name in syutyu.__all__ else f'{module}.{name}'
            for module, name, _value in objects.values()
        }
        listed = set(re.findall(r'^\.\. (?:autoclass|autofunction)::\s+(\S+)$',
                                REFERENCE.read_text(encoding='utf-8'), flags=re.M))
        self.assertEqual(listed, expected)

    def test_no_custom_section_labels(self):
        objects, members = collect()
        for name, value in {**{key: item[2] for key, item in objects.items()}, **members}.items():
            with self.subTest(name=name):
                for line in documented(value).splitlines():
                    self.assertNotIn(line.strip(), CUSTOM_LABELS)

    def test_parameters_are_in_args(self):
        objects, members = collect()
        for name, value in {**{key: item[2] for key, item in objects.items()}, **members}.items():
            if not inspect.isfunction(value):
                continue
            parameters = [item for item in inspect.signature(value).parameters.values()
                          if item.name != 'self']
            if not parameters:
                continue
            doc = documented(value).splitlines()
            with self.subTest(name=name):
                self.assertIn('Args:', [line.strip() for line in doc])
                body = section_lines(doc, 'Args:')
                for parameter in parameters:
                    if parameter.kind is inspect.Parameter.VAR_KEYWORD:
                        prefix = '**'
                    elif parameter.kind is inspect.Parameter.VAR_POSITIONAL:
                        prefix = '*'
                    else:
                        prefix = ''
                    entry = re.compile(r'^\s+' + re.escape(prefix + parameter.name) + r'(\s\(.*\))?:\s')
                    self.assertTrue(any(entry.match(line) for line in body), parameter.name)

    def test_dataclass_fields_are_in_attributes(self):
        objects, _members = collect()
        for name, (_module, _short, value) in objects.items():
            if not inspect.isclass(value) or not is_dataclass(value):
                continue
            doc = documented(value).splitlines()
            with self.subTest(name=name):
                self.assertIn('Attributes:', [line.strip() for line in doc])
                body = section_lines(doc, 'Attributes:')
                for field in fields(value):
                    entry = re.compile(r'^\s+' + re.escape(field.name) + r'(\s\(.*\))?:\s')
                    self.assertTrue(any(entry.match(line) for line in body), field.name)

    def test_inline_literals_are_rst(self):
        objects, members = collect()
        for name, value in {**{key: item[2] for key, item in objects.items()}, **members}.items():
            doc = documented(value)
            with self.subTest(name=name):
                self.assertIsNone(LONE_BACKTICK.search(doc))
                self.assertIsNone(AFTER_KANA.search(doc))
                self.assertIsNone(BEFORE_KANA.search(doc))


if __name__ == '__main__':
    unittest.main()
