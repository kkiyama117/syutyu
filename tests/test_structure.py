from pathlib import Path
import tempfile
import unittest

import numpy as np

from syutyu import read_trajectory
from syutyu.analysis.msd import mean_squared_displacement
from syutyu.analysis.rdf import radial_distribution
from syutyu.io.native import read_msd
from syutyu.kernels import pair_rdf
from syutyu.report.figures import make_figure
from support import run_cli, tiny_dump, write_dump


class StructureTests(unittest.TestCase):
    def test_reordered_atoms_and_extra_charge_column(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'dump'
            t = read_trajectory(write_dump(path, 'id type mol q x y z',
                                           [['2 2 2 -1 9.8 0 0', '1 1 1 1 .2 0 0'],
                                            ['1 1 1 1 .2 0 0', '2 2 2 -1 9.8 0 0']]))
            np.testing.assert_array_equal(t.ids, [1, 2])
            np.testing.assert_array_equal(t.select('type:1'), [True, False])
            np.testing.assert_array_equal(t.select('mol:2'), [False, True])
            with self.assertRaisesRegex(ValueError, 'no atoms'):
                t.select('type:99')
            with self.assertRaisesRegex(ValueError, 'Unavailable'):
                t.select('element:Na')

    def test_fluctuating_volume_rdf_normalization_and_cn(self):
        with tempfile.TemporaryDirectory() as tmp:
            t = read_trajectory(write_dump(
                Path(tmp) / 'dump', 'id type x y z',
                [['1 1 .2 0 0', '2 2 9.8 0 0'], ['1 1 .24 0 0', '2 2 11.76 0 0']],
                boxes=[((0, 10),) * 3, ((0, 12),) * 3]))
            result = radial_distribution(t, [('type:1', 'type:2')], rmax=1, bins=2,
                                         normalization='population')
            shell = 4*np.pi/3*.5**3
            self.assertAlmostEqual(result.tables['rdf']['g_0'][0], (1000+1728)/2/shell)
            np.testing.assert_allclose(result.tables['rdf']['CN_0'], [1, 1])
            with self.assertRaisesRegex(ValueError, 'Variable-cell'):
                mean_squared_displacement(t, ['type:1'], dt_ps=.001)

    def test_same_selection_matches_original_rdf(self):
        with tempfile.TemporaryDirectory() as tmp:
            dump = tiny_dump(Path(tmp))
            t = read_trajectory(dump)
            result = radial_distribution(t, [('all', 'all')], rmax=4, bins=8,
                                         normalization='population')
            expected, cn, _ = pair_rdf(t.positions, t.boxes[0], np.linspace(0, 4, 9),
                                       normalization='population')
            np.testing.assert_allclose(result.tables['rdf']['g_0'], expected)
            np.testing.assert_allclose(result.tables['rdf']['CN_0'], cn)

    def test_intermolecular_exclusion(self):
        with tempfile.TemporaryDirectory() as tmp:
            t = read_trajectory(write_dump(Path(tmp) / 'dump', 'id type mol x y z',
                                           [['1 1 1 1 0 0', '2 2 1 2 0 0']]))
            result = radial_distribution(t, [('type:1', 'type:2')], rmax=3, bins=3,
                                         intermolecular=True, normalization='population')
            self.assertTrue((result.tables['rdf']['CN_0'] == 0).all())
            t.labels['mol'][:] = '0'
            result = radial_distribution(t, [('type:1', 'type:2')], rmax=3, bins=3,
                                         intermolecular=True, normalization='population')
            self.assertEqual(result.tables['rdf']['CN_0'][-1], 1)

    def test_scaled_coordinates_and_image_flags(self):
        with tempfile.TemporaryDirectory() as tmp:
            t = read_trajectory(write_dump(
                Path(tmp) / 'dump', 'id type xs ys zs ix iy iz',
                [['1 1 .9 0 0 0 0 0'], ['1 1 .1 0 0 1 0 0'], ['1 1 .3 0 0 1 0 0']]))
            result = mean_squared_displacement(t, ['all'], dt_ps=.1)
            np.testing.assert_allclose(result.tables['msd']['lag_ps'], [0, 1, 2])
            np.testing.assert_allclose(result.tables['msd']['MSD_0_A2'], [0, 4, 16])

    def test_npt_explicit_unwrapped_requires_lab_opt_in(self):
        with tempfile.TemporaryDirectory() as tmp:
            t = read_trajectory(write_dump(Path(tmp) / 'dump', 'id type xu yu zu',
                                           [['1 1 9 0 0'], ['1 1 12 0 0']],
                                           boxes=[((0, 10),) * 3, ((0, 12),) * 3]))
            with self.assertRaisesRegex(ValueError, 'laboratory'):
                mean_squared_displacement(t, ['type:1'], dt_ps=.1)
            result = mean_squared_displacement(t, ['type:1'], dt_ps=.1, laboratory=True)
            np.testing.assert_allclose(result.tables['msd']['MSD_0_A2'], [0, 9])
            self.assertIn('laboratory', result.summary['coordinate_frame'])

    def test_explicit_unwrapped_preferred_over_wrapped(self):
        with tempfile.TemporaryDirectory() as tmp:
            t = read_trajectory(write_dump(Path(tmp) / 'dump', 'id type x y z xu yu zu',
                                           [['1 1 1 0 0 1 0 0'], ['1 1 7 0 0 17 0 0']]))
            result = mean_squared_displacement(t, ['all'], dt_ps=.1)
            self.assertEqual(result.tables['msd']['MSD_0_A2'][1], 256)

    def test_native_msd_tables(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'Na.msd'
            path.write_text('# TimeStep msd\n0 0\n5 2\n10 4\n')
            result = read_msd({'Na': path}, dt_ps=.001, origin_step=0)
            np.testing.assert_allclose(result.tables['msd']['lag_ps'], [0, .005, .01])
            np.testing.assert_allclose(result.tables['msd']['MSD_0_A2'], [0, 2, 4])
            self.assertEqual(make_figure(result).axes[0].lines[0].get_label(), 'Na')
            with self.assertRaises(ValueError):
                read_msd({'Na': path}, dt_ps=.001, origin_step=0, column=0)

    def test_cli_partial_rdf_and_plot(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            dump = tiny_dump(tmp)
            process = run_cli('rdf', dump, '--pair', 'element:Ar', 'all',
                              '--rmax', 4, '--bins', 8, '--output', tmp / 'result')
            self.assertEqual(process.returncode, 0, process.stderr)
            process = run_cli('plot', tmp / 'result', '--view', 'cn', '--output', tmp / 'cn.png')
            self.assertEqual(process.returncode, 0, process.stderr)
            self.assertTrue((tmp/'cn.png').is_file())


if __name__ == '__main__':
    unittest.main()
