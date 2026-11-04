from pathlib import Path
import tempfile
import unittest

from syutyu import Site, write_charge_data

INPUTS = {
    'lengths': (4.115, 4.115, 4.115),
    'sites': [Site(1, 1.0, (0.0, 0.0, 0.0)), Site(2, -1.0, (0.5, 0.5, 0.5))],
    'masses': {1: 132.90543, 2: 35.4527},
}


class ChargeDataTests(unittest.TestCase):
    def test_two_site_cell_text_is_exact(self):
        expected = (
            'CsCl unit cell\n'
            '\n'
            '2 atoms\n'
            '2 atom types\n'
            '\n'
            '0.0 4.1150000000 xlo xhi\n'
            '0.0 4.1150000000 ylo yhi\n'
            '0.0 4.1150000000 zlo zhi\n'
            '\n'
            'Masses\n'
            '\n'
            '1 132.905430\n'
            '2 35.452700\n'
            '\n'
            'Atoms # charge\n'
            '\n'
            '1 1 +1.000000 0.0000000000 0.0000000000 0.0000000000\n'
            '2 2 -1.000000 2.0575000000 2.0575000000 2.0575000000\n'
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = write_charge_data(Path(tmp) / 'cell.data', title='CsCl unit cell',
                                     **INPUTS)
            self.assertEqual(path.read_text(encoding='utf-8'), expected)

    def test_types_listed_in_ascending_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = write_charge_data(
                Path(tmp) / 'cell.data', (2.0, 3.0, 4.0),
                [Site(2, -1.0, (0.5, 0.5, 0.5)), Site(1, 1.0, (0.0, 0.0, 0.0))],
                {2: 35.4527, 1: 132.90543})
            lines = path.read_text(encoding='utf-8').splitlines()
            masses = lines[lines.index('Masses') + 2:lines.index('Atoms # charge') - 1]
            self.assertEqual([line.split()[0] for line in masses], ['1', '2'])
            atoms = lines[lines.index('Atoms # charge') + 2:]
            self.assertEqual([line.split()[:2] for line in atoms],
                             [['1', '2'], ['2', '1']])

    def test_rejects_bad_inputs(self):
        cases = [
            ('zero length', {'lengths': (0.0, 1.0, 1.0)}),
            ('infinite length', {'lengths': (1.0, float('inf'), 1.0)}),
            ('two lengths', {'lengths': (1.0, 1.0)}),
            ('empty sites', {'sites': []}),
            ('zero type', {'sites': [Site(0, 1.0, (0.0, 0.0, 0.0))]}),
            ('outside cell', {'sites': [Site(1, 1.0, (1.0, 0.0, 0.0))]}),
            ('two coordinates', {'sites': [Site(1, 1.0, (0.0, 0.0))]}),
            ('nan charge', {'sites': [Site(1, float('nan'), (0.0, 0.0, 0.0))]}),
            ('zero mass', {'masses': {1: 0.0, 2: 35.4527}}),
            ('type without mass', {'masses': {1: 132.90543}}),
            ('unused mass', {'masses': {1: 132.90543, 2: 35.4527, 3: 1.0}}),
            ('two-line title', {'title': 'a\nb'}),
        ]
        with tempfile.TemporaryDirectory() as tmp:
            for name, changes in cases:
                with self.subTest(name):
                    with self.assertRaises(ValueError):
                        write_charge_data(Path(tmp) / 'cell.data', **(INPUTS | changes))
                    self.assertEqual(list(Path(tmp).iterdir()), [])

    def test_existing_file_needs_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'cell.data'
            write_charge_data(path, title='first', **INPUTS)
            with self.assertRaises(FileExistsError):
                write_charge_data(path, title='second', **INPUTS)
            write_charge_data(path, title='second', overwrite=True, **INPUTS)
            self.assertEqual(path.read_text(encoding='utf-8').splitlines()[0], 'second')

    def test_exported_from_package(self):
        from syutyu import Site, write_charge_data
        self.assertEqual(Site(1, 1.0, (0.0, 0.0, 0.0)).fractional, (0.0, 0.0, 0.0))
        self.assertTrue(callable(write_charge_data))


if __name__ == '__main__':
    unittest.main()
