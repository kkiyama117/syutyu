"""Feature-package contracts for the lecture's real LAMMPS output formats."""
import dataclasses
from pathlib import Path
import tempfile
import unittest

import numpy as np

from syutyu import (DiffusionFit, fit_diffusion, load_result, read_thermo, save_figure,
                    save_result, summarize_thermo)
from syutyu.io.native import read_msd, read_rdf


class PostprocessingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def rdf_file(self, content):
        path = self.root / 'Ar.gr'
        path.write_text(content)
        return path

    def test_vector_rdf_average_selection_and_result_roundtrip(self):
        path = self.rdf_file('# RDF blocks\n400 2\n1 0.025 0 0\n2 0.075 2 1\n'
                             '500 2\n1 0.025 0 0\n2 0.075 4 3\n')
        result = read_rdf(path, pairs=[('Ar', 'Ar')])
        table = result.tables['rdf']
        np.testing.assert_allclose(table['r_upper_A'].to_numpy(), [.05, .10])
        np.testing.assert_allclose(table['g_0'].to_numpy(), [0, 3])
        np.testing.assert_allclose(table['CN_0'].to_numpy(), [0, 2])
        self.assertEqual(result.summary['block_steps'], [400, 500])
        selected = read_rdf(path, start_step=500)
        np.testing.assert_allclose(selected.tables['rdf']['g_0'].to_numpy(), [0, 4])
        save_result(result, self.root / 'result')
        loaded = load_result(self.root / 'result')
        np.testing.assert_allclose(loaded.tables['rdf']['g_0'].to_numpy(), [0, 3])
        save_figure(result, self.root / 'rdf.png')
        save_figure(result, self.root / 'cn.png', view='cn')
        self.assertGreater((self.root / 'rdf.png').stat().st_size, 1000)

    def test_partial_rdf_pairs_are_not_inferred(self):
        path = self.rdf_file('100 1\n1 0.5 1 2 3 4\n')
        result = read_rdf(path, pairs=[('Na', 'Na'), ('Na', 'Cl')])
        self.assertEqual(result.summary['rdf_pairs'][1]['neighbor'], 'Cl')
        np.testing.assert_allclose(result.tables['rdf']['g_1'].to_numpy(), [3])
        with self.assertRaises(ValueError):
            read_rdf(path, pairs=[('Ar', 'Ar')])

    def test_malformed_and_truncated_rdf_rejected(self):
        cases = ['', '100 2\n1 .025 0 0\n',
                 '100 1\n2 .025 0 0\n', '100 1\n1 .025 nan 0\n',
                 '100 1\n1 .025 0 0\n100 1\n1 .025 0 0\n',
                 '100 2\n1 .025 0 0\n2 .1 1 2\n',
                 '100 1\n1 .025 0 0\n200 1\n1 .05 0 0\n']
        for text in cases:
            with self.subTest(text=text), self.assertRaises(ValueError):
                read_rdf(self.rdf_file(text))
        with self.assertRaises(ValueError):
            read_rdf(self.rdf_file('100 1\n1 .025 0 0\n'), start_step=200)

    def test_msd_origin_is_compute_step_not_first_saved_step(self):
        path = self.root / 'Ar.msd'
        path.write_text('# LAMMPS compute msd\n400 .4 .4 .4 1.2\n500 .8 .8 .8 2.4\n600 1.2 1.2 1.2 3.6\n')
        result = read_msd({'Ar': path}, dt_ps=.01, origin_step=300)
        table = result.tables['msd']
        np.testing.assert_allclose(table['lag_ps'].to_numpy(), [1, 2, 3])
        fit = fit_diffusion(table['lag_ps'].to_numpy(), table['MSD_0_A2'].to_numpy(), window_ps=(1, 3))
        self.assertAlmostEqual(fit.D_m2_s, 2e-9)
        with self.assertRaises(TypeError):
            read_msd({'Ar': path}, dt_ps=.01)

    def test_diffusion_intercept_units_and_explicit_window(self):
        time = np.arange(31.)
        msd = 6 * .243 * time + .2  # Rahman printed asymptote, synthetic test only
        msd[:5] = .1 * time[:5] ** 2
        fit = fit_diffusion(time, msd, window_ps=(10, 30))
        self.assertAlmostEqual(fit.D_m2_s, 2.43e-9, places=18)
        self.assertAlmostEqual(fit.regression.intercept, .2)
        self.assertEqual(fit.samples, 21)
        np.testing.assert_allclose(fit.regression.intercept + fit.regression.slope * np.array([10, 30]),
                                   msd[[10, 30]])
        for window in ((-1, 30), (10, 31), (20, 10), (10, 11)):
            with self.subTest(window=window), self.assertRaises(ValueError):
                fit_diffusion(time, msd, window_ps=window)
        with self.assertRaises(ValueError):
            fit_diffusion([0, 1, 2], [2, 1, 0], window_ps=(0, 2))
        with self.assertRaises(ValueError):
            fit_diffusion([0, 1, 2], [0, np.nan, 2], window_ps=(0, 2))

    def test_diffusion_fit_holds_linregress(self):
        fit = fit_diffusion(np.arange(5.), 6 * np.arange(5.) + 1, window_ps=(0, 4))
        self.assertAlmostEqual(fit.regression.slope, 6)
        self.assertAlmostEqual(fit.regression.intercept, 1)
        self.assertAlmostEqual(fit.regression.rvalue ** 2, 1)
        self.assertAlmostEqual(fit.D_m2_s, 1e-8, places=20)
        self.assertEqual([field.name for field in dataclasses.fields(DiffusionFit)],
                         ['window_ps', 'samples', 'regression', 'D_m2_s', 'limitation'])
        self.assertFalse(hasattr(fit, 'predict'))

    def test_separate_thermo_blocks_and_density_fluctuations(self):
        path = self.root / 'log.lammps'
        path.write_text('Step Temp Density\n0 500 1.374\nLoop time of 1\n'
                        'Step Temp Density\n300 94 1.3\n400 96 1.4\n500 95 1.5\nLoop time of 1\n')
        table = read_thermo(path, block=-1, dt_ps=.01)
        stats = summarize_thermo(table, start_ps=4, stop_ps=5)
        rows = {row['column']: row for row in stats.iter_rows(named=True)}
        self.assertEqual(rows['Temp']['mean'], 95.5)
        self.assertAlmostEqual(rows['Density']['mean'], 1.45)
        self.assertEqual(rows['Density']['samples'], 2)
        self.assertAlmostEqual(rows['Density']['std'], .05)
        with self.assertRaises(ValueError):
            summarize_thermo(table, start_ps=10)


if __name__ == '__main__':
    unittest.main()
