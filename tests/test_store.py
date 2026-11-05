"""schema3 の表が polars の DataFrame になり、保存して読み戻せること。"""
from pathlib import Path
import tempfile
import unittest

import numpy as np
import polars as pl

from syutyu import Result, load_result, msd, save_result
from syutyu.io.native import read_msd
from support import memory_universe


def sliding_trajectory():
    positions = np.array([[[0., 0, 0], [5., 0, 0]], [[1., 0, 0], [5., 0, 0]],
                          [[2., 0, 0], [5., 0, 0]]])
    return memory_universe([0, 10, 20], positions, np.full((3, 3), 10.), ['1', '2'])


class ResultStoreTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_result_tables_are_polars(self):
        path = self.root / 'Ar.msd'
        path.write_text('400 .4 .4 .4 1.2\n500 .8 .8 .8 2.4\n')
        table = read_msd({'Ar': path}, dt_ps=.01, origin_step=300).tables['msd']
        self.assertIsInstance(table, pl.DataFrame)
        self.assertEqual(table.columns, ['lag_ps', 'MSD_0_A2'])
        self.assertEqual(table.schema, {'lag_ps': pl.Float64, 'MSD_0_A2': pl.Float64})
        self.assertEqual(table['MSD_0_A2'].to_list(), [1.2, 2.4])

    def test_integer_columns_are_int64(self):
        result = msd(sliding_trajectory(), method='custom', dt_ps=.001, coordinates='unwrapped')
        table = result.tables['msd']
        self.assertEqual(table.schema['origin_atom_pairs_0'], pl.Int64)
        self.assertEqual(table.schema['lag_ps'], pl.Float64)
        self.assertEqual(table.schema['MSD_0_A2'], pl.Float64)

    def test_roundtrip_is_exact(self):
        result = msd(sliding_trajectory(), method='custom', dt_ps=.001, coordinates='unwrapped')
        output = save_result(result, self.root / 'msd-result')
        loaded = load_result(output)
        self.assertTrue(loaded.tables['msd'].equals(result.tables['msd']))
        self.assertEqual(loaded.tables['msd'].schema, result.tables['msd'].schema)
        self.assertEqual(loaded.summary, {**result.summary, 'schema_version': 3, 'tables': ['msd']})

    def test_non_polars_table_rejected(self):
        output = self.root / 'plain-dict'
        with self.assertRaises(ValueError) as caught:
            save_result(Result({'msd': {'lag_ps': np.array([0.])}}, {}), output)
        self.assertEqual(str(caught.exception), 'Result tables must be polars DataFrames')
        self.assertFalse(output.exists())

    def test_numpy_written_csv_still_loads(self):
        directory = self.root / 'numpy-result'
        directory.mkdir()
        (directory / 'msd.csv').write_text(
            'lag_ps,MSD_0_A2\n0.000000000000000000e+00,1.000000000000000000e+00\n'
            '1.000000000000000000e+00,2.500000000000000000e+00\n')
        (directory / 'summary.json').write_text('{"schema_version": 3, "tables": ["msd"]}')
        table = load_result(directory).tables['msd']
        self.assertEqual(table['lag_ps'].to_list(), [0.0, 1.0])
        self.assertEqual(table['MSD_0_A2'].to_list(), [1.0, 2.5])
        self.assertEqual(dict(table.schema), {'lag_ps': pl.Float64, 'MSD_0_A2': pl.Float64})

    def test_empty_csv_rejected(self):
        directory = self.root / 'header-only'
        directory.mkdir()
        (directory / 'msd.csv').write_text('lag_ps,MSD_0_A2\n')
        (directory / 'summary.json').write_text('{"schema_version": 3, "tables": ["msd"]}')
        with self.assertRaises(ValueError) as caught:
            load_result(directory)
        self.assertTrue(str(caught.exception).startswith('Empty or unnamed CSV'))

    def test_zero_row_table_refused(self):
        output = self.root / 'empty-table'
        table = pl.DataFrame({'lag_ps': []}, schema={'lag_ps': pl.Float64})
        with self.assertRaises(ValueError) as caught:
            save_result(Result({'msd': table}, {}), output)
        self.assertEqual(str(caught.exception), 'Result tables must have rows')
        self.assertFalse(output.exists())


if __name__ == '__main__':
    unittest.main()
