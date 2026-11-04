"""Method dispatch, finite-N RDF, native VACF and integrated output contracts."""
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from syutyu import Trajectory, load_result, msd, rdf, save_result, vacf
from syutyu.analysis.rdf import radial_distribution
from syutyu.io.native import read_vacf
from syutyu.kernels import pair_rdf
from syutyu.report.figures import available_views, make_figure
from support import run_cli, write_dump


def tiny_trajectory():
    positions = np.array([[[0., 0, 0], [2., 0, 0]], [[1., 0, 0], [3., 0, 0]],
                          [[3., 0, 0], [5., 0, 0]]])
    velocities = np.array([[[1., 0, 0], [2., 0, 0]], [[-1., 0, 0], [1., 0, 0]],
                           [[2., 0, 0], [-2., 0, 0]]])
    return Trajectory(steps=np.array([0, 10, 20]), ids=np.array([1, 2]), positions=positions,
                      boxes=np.full((3, 3), 10.), origins=np.zeros((3, 3)),
                      labels={'type': np.array(['1', '1'])}, unwrapped=positions.copy(),
                      unwrap_source='xu/yu/zu', velocities=velocities)


class MethodChoiceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_custom_finite_rdf_matches_colab_unique_pair_formula(self):
        t = tiny_trajectory()
        result = rdf(t, method='custom', pairs=[('all', 'all')], rmax=4, bins=4)
        expected = np.zeros(4)
        expected[2] = 1000 / (4 * np.pi / 3 * (3**3 - 2**3))  # H=1; N(N-1)/2=1
        np.testing.assert_allclose(result.tables['rdf']['g_0'], expected)
        population = radial_distribution(t, [('all', 'all')], rmax=4, bins=4,
                                         normalization='population')
        np.testing.assert_allclose(expected, 2 * population.tables['rdf']['g_0'])
        kernel = pair_rdf(t.positions, t.boxes[0], np.arange(5.), normalization='finite')[0]
        np.testing.assert_allclose(kernel, expected)
        self.assertEqual(result.summary['rdf_normalization'], 'finite')
        self.assertEqual(result.summary['method'], 'custom')
        np.testing.assert_allclose(result.tables['rdf']['CN_0'], [0, 0, 1, 1])

    def test_last_radial_edge_is_excluded_as_in_colab_and_lammps(self):
        t = tiny_trajectory()
        kernel, cn, counts = pair_rdf(t.positions, t.boxes[0], np.array([0., 1., 2.]), normalization='finite')
        np.testing.assert_array_equal(counts, [0, 0])  # all pair distances equal rmax=2
        result = rdf(t, method='custom', pairs=[('all', 'all')], rmax=2, bins=2)
        np.testing.assert_array_equal(result.tables['rdf']['g_0'], kernel)
        np.testing.assert_array_equal(result.tables['rdf']['CN_0'], cn)

    def test_native_rdf_dispatch_does_not_reconstruct_coordinates(self):
        path = self.root / 'native.gr'
        path.write_text('100 2\n1 .5 0 0\n2 1.5 2 1\n')
        result = rdf(path, method='lammps')
        np.testing.assert_allclose(result.tables['rdf']['g_0'], [0, 2])
        self.assertEqual(result.summary['rdf_normalization'], 'finite')
        with self.assertRaises(ValueError):
            rdf(path, method='lammps', normalization='population')
        with self.assertRaises(ValueError):
            rdf(path, method='lammps', rmax=2)
        with self.assertRaises(ValueError):
            rdf(path, method='automatic')

    def test_custom_msd_single_and_multiple_origins_are_distinct(self):
        t = tiny_trajectory()
        single = msd(t, method='custom', dt_ps=.1, origin_mode='single')
        multiple = msd(t, method='custom', dt_ps=.1)
        np.testing.assert_allclose(single.tables['msd']['MSD_0_A2'], [0, 1, 9])
        np.testing.assert_allclose(multiple.tables['msd']['MSD_0_A2'], [0, 2.5, 9])
        self.assertEqual(single.summary['origin_mode'], 'single')
        self.assertEqual(multiple.summary['origin_mode'], 'multiple')
        with self.assertRaises(ValueError):
            msd(t, method='custom', dt_ps=.1, origin_step=0)

    def test_native_msd_requires_origin_and_rejects_multiple(self):
        path = self.root / 'native.msd'
        path.write_text('100 0 0 0 0\n110 1 0 0 1\n120 4 0 0 4\n')
        result = msd(path, method='lammps', dt_ps=.01, origin_step=100)
        np.testing.assert_allclose(result.tables['msd']['lag_ps'], [0, .1, .2])
        self.assertEqual(result.summary['origin_mode'], 'single')
        with self.assertRaises(ValueError):
            msd(path, method='lammps', dt_ps=.01)
        with self.assertRaises(ValueError):
            msd(path, method='lammps', dt_ps=.01, origin_step=100, origin_mode='multiple')

    def vacf_file(self):
        path = self.root / 'native.vacf'
        path.write_text('# timestep xx yy zz sum\n100 1 0 0 1\n110 -1 0 0 -1\n120 .5 .5 0 1\n')
        return path

    def test_native_vacf_units_negative_values_and_normalization(self):
        result = vacf(self.vacf_file(), method='lammps', dt_ps=.01,
                                origin_step=100, velocity_unit='A/fs', max_lag_ps=.1)
        np.testing.assert_allclose(result.tables['vacf']['VACF_0_A2_ps2'], [1e6, -1e6])
        np.testing.assert_allclose(result.tables['vacf']['VACF_0_normalized'], [1, -1])
        np.testing.assert_allclose(result.tables['vacf']['lag_ps'], [0, .1])
        self.assertEqual(result.summary['origin_mode'], 'single')
        self.assertEqual(result.summary['method'], 'lammps')
        with self.assertRaises(ValueError):
            vacf(self.vacf_file(), method='lammps', dt_ps=.01,
                           origin_step=100, velocity_unit='A/fs', origin_mode='multiple')

    def test_native_vacf_rejects_missing_origin_bad_components_and_wrong_units(self):
        for content in ('110 1 0 0 1\n120 -1 0 0 -1\n',
                        '100 1 0 0 2\n110 -1 0 0 -1\n',
                        '100 0 0 0 0\n110 -1 0 0 -1\n',
                        '100 nan 0 0 nan\n110 -1 0 0 -1\n'):
            path = self.root / 'bad.vacf'
            path.write_text(content)
            with self.subTest(content=content), self.assertRaises(ValueError):
                read_vacf(path, dt_ps=.01, origin_step=100, velocity_unit='A/fs')
        with self.assertRaises(ValueError):
            read_vacf(self.vacf_file(), dt_ps=.01, origin_step=100, velocity_unit='m/s')

    def test_vacf_save_load_and_both_plot_views(self):
        result = vacf(tiny_trajectory(), method='custom', dt_ps=.1, velocity_unit='A/ps')
        save_result(result, self.root / 'result')
        loaded = load_result(self.root / 'result')
        self.assertEqual(set(available_views(loaded)), {'vacf', 'vacf-normalized'})
        raw = make_figure(loaded, 'vacf').axes[0].lines[0].get_ydata()
        normalized = make_figure(loaded, 'vacf-normalized').axes[0].lines[0].get_ydata()
        np.testing.assert_allclose(normalized, raw / raw[0])
        self.assertEqual(loaded.summary['origin_mode'], 'multiple')

    def test_cli_native_and_custom_mode_errors_are_explicit(self):
        output = self.root / 'cli-native'
        native = ('vacf', self.vacf_file(), '--method', 'lammps', '--dt-ps', '.01',
                  '--origin-step', '100', '--velocity-unit', 'A/ps')
        proc = run_cli(*native, '--output', output)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(load_result(output).summary['method'], 'lammps')
        proc = run_cli(*native, '--origins', 'multiple', '--output', self.root / 'bad')
        self.assertEqual(proc.returncode, 2)
        self.assertIn('single origin', proc.stderr)


if __name__ == '__main__':
    unittest.main()
