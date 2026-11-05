from pathlib import Path
import math
import tempfile
import unittest

import numpy as np
import polars as pl

from syutyu import read_thermo, read_thermo_blocks, summarize_thermo


class ThermoTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.path = Path(self.temporary.name) / 'log.lammps'

    def log(self, text):
        self.path.write_text(text)
        return self.path

    def test_multiple_blocks_different_headers_and_reset_steps(self):
        path = self.log('''LAMMPS synthetic fixture
Step Temp PotEng
100 300 -2.5
WARNING: sample warning
200 310 -2.4
Loop time of 0.1 on 1 procs
reset_timestep 0
Step Press Volume
0 1.0 1234
5 1.1 1235
Loop time of 0.2 on 1 procs
Total wall time: 0:00:00
''')
        blocks = read_thermo_blocks(path)
        self.assertEqual(len(blocks), 2)
        self.assertEqual(blocks[0].columns, ['Step', 'Temp', 'PotEng'])
        self.assertEqual(blocks[1].columns, ['Step', 'Press', 'Volume'])
        np.testing.assert_array_equal(blocks[0]['Step'].to_numpy(), [100, 200])
        np.testing.assert_array_equal(blocks[1]['Step'].to_numpy(), [0, 5])
        table = read_thermo(path, dt_ps=.002)
        np.testing.assert_allclose(table['time_ps'].to_numpy(), [0, .01])
        self.assertNotIn('Temp', table.columns)
        first = read_thermo(path, block=0, dt_ps=.001)
        np.testing.assert_allclose(first['time_ps'].to_numpy(), [.1, .2])
        np.testing.assert_array_equal(first['PotEng'].to_numpy(), [-2.5, -2.4])
        np.testing.assert_array_equal(read_thermo(path, block=-2, dt_ps=1)['Step'].to_numpy(), [100, 200])

    def test_blocks_are_polars_tables(self):
        path = self.log('''LAMMPS synthetic fixture
Step Temp PotEng
100 300 -2.5
WARNING: sample warning
200 310 -2.4
Loop time of 0.1 on 1 procs
reset_timestep 0
Step Press Volume
0 1.0 1234
5 1.1 1235
Loop time of 0.2 on 1 procs
Total wall time: 0:00:00
''')
        blocks = read_thermo_blocks(path)
        self.assertEqual(len(blocks), 2)
        for index, names in enumerate((['Step', 'Temp', 'PotEng'], ['Step', 'Press', 'Volume'])):
            with self.subTest(block=index):
                self.assertIsInstance(blocks[index], pl.DataFrame)
                self.assertEqual(blocks[index].columns, names)
                self.assertEqual(blocks[index].dtypes, [pl.Float64] * len(names))
        self.assertEqual(read_thermo(path, dt_ps=.01).columns[0], 'time_ps')

    def test_summary_is_a_polars_table(self):
        table = pl.DataFrame({'time_ps': [0., 1., 2.], 'Step': [0, 10, 20],
                              'Temp': [90., 100., 110.], 'Press': [1., 1., 1.]})
        stats = summarize_thermo(table)
        self.assertEqual(stats.columns, ['column', 'mean', 'std', 'samples'])
        self.assertEqual(stats.dtypes, [pl.String, pl.Float64, pl.Float64, pl.Int64])
        self.assertEqual(stats['column'].to_list(), ['Temp', 'Press'])
        rows = {row['column']: row for row in stats.iter_rows(named=True)}
        self.assertEqual(rows['Temp']['mean'], 100.0)
        self.assertAlmostEqual(rows['Temp']['std'], math.sqrt(200 / 3))
        self.assertEqual(rows['Temp']['samples'], 3)

    def test_summary_rejects_a_dict(self):
        with self.assertRaises(TypeError) as raised:
            summarize_thermo({'time_ps': np.array([0.]), 'Temp': np.array([1.])})
        self.assertEqual(str(raised.exception), 'table must be a polars DataFrame')

    def test_identical_headers_do_not_merge_runs(self):
        path = self.log('Step Temp\n0 100\n1 101\nLoop time of 0\n'
                        'Step Temp\n0 200\n1 201\nLoop time of 0\n')
        self.assertEqual(len(read_thermo_blocks(path)), 2)
        np.testing.assert_array_equal(read_thermo(path, dt_ps=1)['Temp'].to_numpy(), [200, 201])

    def test_single_row_block_is_allowed(self):
        path = self.log('Step Temp\n17 298\nLoop time of 0\n')
        np.testing.assert_array_equal(read_thermo(path, dt_ps=2)['time_ps'].to_numpy(), [34])

    def test_dt_ps_is_required_finite_and_positive(self):
        path = self.log('Step Temp\n0 298\nLoop time of 0\n')
        with self.assertRaises(TypeError):
            read_thermo(path)
        for dt in (0, -1, np.nan, np.inf, -np.inf):
            with self.subTest(dt=dt), self.assertRaisesRegex(ValueError, 'dt_ps'):
                read_thermo(path, dt_ps=dt)

    def test_block_index_out_of_range(self):
        path = self.log('Step Temp\n0 298\nLoop time of 0\n')
        for index in (1, -2, 100):
            with self.subTest(index=index), self.assertRaises(IndexError):
                read_thermo(path, block=index, dt_ps=1)
        with self.assertRaises(ValueError):
            read_thermo(path, block=.5, dt_ps=1)

    def test_finite_validation_only_applies_to_selected_block(self):
        path = self.log('Step Temp\n0 nan\n1 inf\nLoop time of 0\n'
                        'Step Temp\n0 298\n1 299\nLoop time of 0\n')
        self.assertEqual(len(read_thermo_blocks(path)), 2)
        read_thermo(path, dt_ps=1)
        with self.assertRaisesRegex(ValueError, 'non-finite'):
            read_thermo(path, block=0, dt_ps=1)
        path = self.log('Step Temp\nnan 298\nLoop time of 0\n')
        with self.assertRaisesRegex(ValueError, 'non-finite'):
            read_thermo(path, dt_ps=1)

    def test_duplicate_and_decreasing_steps_rejected(self):
        for steps in ('0 298\n0 299', '1 298\n0 299'):
            with self.subTest(steps=steps):
                path = self.log(f'Step Temp\n{steps}\nLoop time of 0\n')
                with self.assertRaisesRegex(ValueError, 'ascending'):
                    read_thermo(path, dt_ps=1)

    def test_empty_and_absent_blocks(self):
        path = self.log('LAMMPS with no run\n')
        self.assertEqual(read_thermo_blocks(path), [])
        with self.assertRaisesRegex(ValueError, 'No thermo blocks'):
            read_thermo(path, dt_ps=1)
        path = self.log('Step Temp\nLoop time of 0\n')
        self.assertEqual(read_thermo_blocks(path)[0]['Step'].to_numpy().shape, (0,))
        with self.assertRaisesRegex(ValueError, 'empty'):
            read_thermo(path, dt_ps=1)

    def test_malformed_or_unfinished_blocks_rejected(self):
        for text in ('Step Temp\n0\nLoop time of 0\n',
                     'Step Temp\n0 bad\nLoop time of 0\n',
                     'Step Step\n0 0\nLoop time of 0\n',
                     'Step Temp\n0 298\n',
                     'Step Temp\n0 298\nStep Temp\n1 299\nLoop time of 0\n'):
            with self.subTest(text=text), self.assertRaises(ValueError):
                read_thermo_blocks(self.log(text))

    def test_derived_time_must_be_finite(self):
        path = self.log('Step Temp\n1e308 298\nLoop time of 0\n')
        with self.assertRaisesRegex(ValueError, 'Derived time_ps'):
            read_thermo(path, dt_ps=10)


if __name__ == '__main__':
    unittest.main()
