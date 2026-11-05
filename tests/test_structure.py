from pathlib import Path
import tempfile
import unittest

import MDAnalysis
import numpy as np

import syutyu
from syutyu import read_trajectory
from syutyu.analysis.msd import mean_squared_displacement
from syutyu.analysis.rdf import radial_distribution
from syutyu.io.native import read_msd
from syutyu.kernels import pair_rdf
from syutyu.report.figures import make_figure
from support import frames_of, memory_universe, run_cli, tiny_dump, write_dump


class StructureTests(unittest.TestCase):
    def test_reordered_atoms_and_extra_charge_column(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'dump'
            t = read_trajectory(write_dump(path, 'id type mol q x y z',
                                           [['2 2 2 -1 9.8 0 0', '1 1 1 1 .2 0 0'],
                                            ['1 1 1 1 .2 0 0', '2 2 2 -1 9.8 0 0']]))
            np.testing.assert_array_equal(t.atoms.ids, [1, 2])
            np.testing.assert_array_equal(t.select_atoms('type 1').ids, [1])
            np.testing.assert_array_equal(t.select_atoms('resid 2').ids, [2])

    def test_fluctuating_volume_rdf_normalization_and_cn(self):
        with tempfile.TemporaryDirectory() as tmp:
            t = read_trajectory(write_dump(
                Path(tmp) / 'dump', 'id type x y z',
                [['1 1 .2 0 0', '2 2 9.8 0 0'], ['1 1 .24 0 0', '2 2 11.76 0 0']],
                boxes=[((0, 10),) * 3, ((0, 12),) * 3]))
            result = radial_distribution(t, [('type 1', 'type 2')], rmax=1, bins=2,
                                         normalization='population')
            shell = 4*np.pi/3*.5**3
            self.assertAlmostEqual(result.tables['rdf']['g_0'].to_numpy()[0], (1000+1728)/2/shell)
            np.testing.assert_allclose(result.tables['rdf']['CN_0'].to_numpy(), [1, 1])
            with self.assertRaisesRegex(ValueError, 'Variable-cell'):
                mean_squared_displacement(t, ['type 1'], dt_ps=.001, coordinates='wrapped')

    def test_same_selection_matches_original_rdf(self):
        with tempfile.TemporaryDirectory() as tmp:
            dump = tiny_dump(Path(tmp))
            t = read_trajectory(dump)
            result = radial_distribution(t, [('all', 'all')], rmax=4, bins=8,
                                         normalization='population')
            _steps, positions, boxes, _velocities = frames_of(t)
            expected, cn, _ = pair_rdf(positions, boxes[0], np.linspace(0, 4, 9),
                                       normalization='population')
            np.testing.assert_allclose(result.tables['rdf']['g_0'].to_numpy(), expected)
            np.testing.assert_allclose(result.tables['rdf']['CN_0'].to_numpy(), cn)

    def test_intermolecular_exclusion(self):
        with tempfile.TemporaryDirectory() as tmp:
            t = read_trajectory(write_dump(Path(tmp) / 'dump', 'id type mol x y z',
                                           [['1 1 1 1 0 0', '2 2 1 2 0 0']]))
            result = radial_distribution(t, [('type 1', 'type 2')], rmax=3, bins=3,
                                         intermolecular=True, normalization='population')
            self.assertTrue((result.tables['rdf']['CN_0'].to_numpy() == 0).all())
            t.residues.resids = np.zeros(t.residues.n_residues, dtype=int)
            result = radial_distribution(t, [('type 1', 'type 2')], rmax=3, bins=3,
                                         intermolecular=True, normalization='population')
            self.assertEqual(result.tables['rdf']['CN_0'].to_numpy()[-1], 1)

    def test_scaled_coordinates_and_image_flags(self):
        with tempfile.TemporaryDirectory() as tmp:
            t = read_trajectory(write_dump(
                Path(tmp) / 'dump', 'id type xs ys zs ix iy iz',
                [['1 1 .9 0 0 0 0 0'], ['1 1 .1 0 0 1 0 0'], ['1 1 .3 0 0 1 0 0']]))
            result = mean_squared_displacement(t, ['all'], dt_ps=.1, coordinates='wrapped')
            np.testing.assert_allclose(result.tables['msd']['lag_ps'].to_numpy(), [0, 1, 2])
            np.testing.assert_allclose(result.tables['msd']['MSD_0_A2'].to_numpy(), [0, 4, 16])

    def test_npt_explicit_unwrapped_requires_lab_opt_in(self):
        with tempfile.TemporaryDirectory() as tmp:
            t = read_trajectory(write_dump(Path(tmp) / 'dump', 'id type xu yu zu',
                                           [['1 1 9 0 0'], ['1 1 12 0 0']],
                                           boxes=[((0, 10),) * 3, ((0, 12),) * 3]))
            with self.assertRaisesRegex(ValueError, 'laboratory'):
                mean_squared_displacement(t, ['type 1'], dt_ps=.1, coordinates='unwrapped')
            result = mean_squared_displacement(t, ['type 1'], dt_ps=.1, coordinates='unwrapped',
                                               laboratory=True)
            np.testing.assert_allclose(result.tables['msd']['MSD_0_A2'].to_numpy(), [0, 9])
            self.assertIn('laboratory', result.summary['coordinate_frame'])

    def test_explicit_unwrapped_preferred_over_wrapped(self):
        with tempfile.TemporaryDirectory() as tmp:
            t = read_trajectory(write_dump(Path(tmp) / 'dump', 'id type x y z xu yu zu',
                                           [['1 1 1 0 0 1 0 0'], ['1 1 7 0 0 17 0 0']]))
            result = mean_squared_displacement(t, ['all'], dt_ps=.1, coordinates='unwrapped')
            self.assertEqual(result.tables['msd']['MSD_0_A2'].to_numpy()[1], 256)

    def test_com_drift_removed_for_uniform_translation(self):
        with tempfile.TemporaryDirectory() as tmp:
            t = read_trajectory(write_dump(Path(tmp) / 'dump', 'id type xu yu zu',
                                           [['1 1 0 0 0', '2 2 1 0 0'],
                                            ['1 1 1 0 0', '2 2 2 0 0'],
                                            ['1 1 2 0 0', '2 2 3 0 0']]))
            removed = mean_squared_displacement(t, ['all'], dt_ps=.1, coordinates='unwrapped',
                                                com_masses={'1': 1.0, '2': 3.0})
            np.testing.assert_array_equal(removed.tables['msd']['MSD_0_A2'].to_numpy(), [0, 0, 0])
            self.assertEqual(removed.summary['com_max_displacement_A'], 2.0)
            self.assertTrue(removed.summary['com_drift_removed'])
            plain = mean_squared_displacement(t, ['all'], dt_ps=.1, coordinates='unwrapped')
            np.testing.assert_array_equal(plain.tables['msd']['MSD_0_A2'].to_numpy(), [0, 1, 4])
            self.assertFalse(plain.summary['com_drift_removed'])

    def test_com_is_mass_weighted_over_all_atoms(self):
        with tempfile.TemporaryDirectory() as tmp:
            t = read_trajectory(write_dump(Path(tmp) / 'dump', 'id type xu yu zu',
                                           [['1 1 0 0 0', '2 2 5 0 0'],
                                            ['1 1 3 0 0', '2 2 5 0 0']]))
            result = mean_squared_displacement(t, ['type 1', 'type 2'], dt_ps=.1,
                                               coordinates='unwrapped',
                                               com_masses={'1': 1.0, '2': 3.0})
            np.testing.assert_array_equal(result.tables['msd']['MSD_0_A2'].to_numpy(), [0, 2.25**2])
            np.testing.assert_array_equal(result.tables['msd']['MSD_1_A2'].to_numpy(), [0, .75**2])

    def test_com_masses_errors(self):
        typed = [['1 1 0 0 0', '2 2 1 0 0'], ['1 1 1 0 0', '2 2 1 0 0']]
        plain = [['1 0 0 0', '2 1 0 0'], ['1 1 0 0', '2 1 0 0']]
        rows = [
            ('id xu yu zu', plain, {'1': 1.0}, False,
             'com_masses needs a type column in the dump'),
            ('id type xu yu zu', typed, {'1': 1.0}, False,
             "com_masses has no mass for types: ['2']"),
            ('id type xu yu zu', typed, {'1': 1.0, '2': 0}, False,
             'com_masses: mass of type 2 must be finite and > 0'),
            ('id type xu yu zu', typed, {'1': float('nan'), '2': 1.0}, False,
             'com_masses: mass of type 1 must be finite and > 0'),
            ('id type xu yu zu', typed, {'1': 1.0, '2': 1.0}, True,
             'com_masses cannot be combined with laboratory=True'),
        ]
        with tempfile.TemporaryDirectory() as tmp:
            for number, (columns, frames, com_masses, laboratory, message) in enumerate(rows):
                with self.subTest(row=number):
                    t = read_trajectory(write_dump(Path(tmp) / f'dump-{number}', columns, frames))
                    with self.assertRaises(ValueError) as caught:
                        mean_squared_displacement(t, ['all'], dt_ps=.1, coordinates='unwrapped',
                                                  com_masses=com_masses, laboratory=laboratory)
                    self.assertEqual(str(caught.exception), message)

    def test_com_masses_default_is_unchanged(self):
        with tempfile.TemporaryDirectory() as tmp:
            t = read_trajectory(write_dump(Path(tmp) / 'dump', 'id type xu yu zu',
                                           [['1 1 0 0 0', '2 2 1 0 0'],
                                            ['1 1 1 0 0', '2 2 1 0 0']]))
            omitted = mean_squared_displacement(t, ['all'], dt_ps=.1, coordinates='unwrapped')
            explicit = mean_squared_displacement(t, ['all'], dt_ps=.1, coordinates='unwrapped',
                                                 com_masses=None)
            for column in omitted.tables['msd'].columns:
                np.testing.assert_array_equal(explicit.tables['msd'][column].to_numpy(),
                                              omitted.tables['msd'][column].to_numpy())
            self.assertIs(omitted.summary['com_drift_removed'], False)
            self.assertIs(explicit.summary['com_drift_removed'], False)

    def test_native_msd_tables(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'Na.msd'
            path.write_text('# TimeStep msd\n0 0\n5 2\n10 4\n')
            result = read_msd({'Na': path}, dt_ps=.001, origin_step=0)
            np.testing.assert_allclose(result.tables['msd']['lag_ps'].to_numpy(), [0, .005, .01])
            np.testing.assert_allclose(result.tables['msd']['MSD_0_A2'].to_numpy(), [0, 2, 4])
            self.assertEqual(make_figure(result).axes[0].lines[0].get_label(), 'Na')
            with self.assertRaises(ValueError):
                read_msd({'Na': path}, dt_ps=.001, origin_step=0, column=0)

    def test_cli_partial_rdf_and_plot(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            dump = tiny_dump(tmp)
            process = run_cli('rdf', dump, '--pair', 'element Ar', 'all',
                              '--rmax', 4, '--bins', 8, '--output', tmp / 'result')
            self.assertEqual(process.returncode, 0, process.stderr)
            process = run_cli('plot', tmp / 'result', '--view', 'cn', '--output', tmp / 'cn.png')
            self.assertEqual(process.returncode, 0, process.stderr)
            self.assertTrue((tmp/'cn.png').is_file())

    def test_read_trajectory_returns_universe(self):
        with tempfile.TemporaryDirectory() as tmp:
            t = read_trajectory(write_dump(Path(tmp) / 'dump', 'id type element mol x y z',
                                           [['2 2 Na 2 9 0 0', '1 1 Cl 1 1 0 0']]))
            self.assertIsInstance(t, MDAnalysis.Universe)
            np.testing.assert_array_equal(t.atoms.ids, [1, 2])
            self.assertEqual(list(t.atoms.elements), ['Cl', 'Na'])
            np.testing.assert_array_equal(t.atoms.resids, [1, 2])
            self.assertEqual(list(t.atoms.types), ['1', '2'])

    def test_read_trajectory_rejects_what_mdanalysis_accepts(self):
        triclinic = ('ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\n'
                     'ITEM: BOX BOUNDS xy xz yz pp pp pp\n0 10 1\n0 10 0\n0 10 0\n'
                     'ITEM: ATOMS id type x y z\n1 1 1 1 1\n')
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / 'triclinic').write_text(triclinic)
            rows = [
                (tmp / 'triclinic',
                 'Only orthogonal periodic pp pp pp cells are supported (no triclinic tilt)'),
                (write_dump(tmp / 'vx-only', 'id type x y z vx', [['1 1 0 0 0 1'], ['1 1 1 0 0 1']]),
                 'Velocity output must include all of vx/vy/vz or none of them'),
                (write_dump(tmp / 'nan', 'id type x y z', [['1 1 nan 0 0'], ['1 1 1 0 0']]),
                 'Non-finite coordinates'),
                (write_dump(tmp / 'steps', 'id type x y z', [['1 1 0 0 0'], ['1 1 1 0 0']],
                            steps=[10, 0]),
                 'Need increasing timesteps; split trajectories at reset_timestep'),
            ]
            for path, message in rows:
                with self.subTest(message=message):
                    with self.assertRaises(ValueError) as caught:
                        read_trajectory(path)
                    self.assertEqual(str(caught.exception), message)

    def test_unwrapped_positions_when_the_dump_allows(self):
        with tempfile.TemporaryDirectory() as tmp:
            t = read_trajectory(write_dump(Path(tmp) / 'dump', 'id type x y z ix iy iz',
                                           [['1 1 9 0 0 0 0 0'], ['1 1 1 0 0 1 0 0']]))
            t.trajectory[1]
            self.assertAlmostEqual(float(t.atoms.positions[0][0]), 11.0)

    def test_msd_coordinates_is_required_and_checked(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = write_dump(Path(tmp) / 'dump', 'id type x y z',
                              [['1 1 0 0 0'], ['1 1 1 0 0']])
            with self.assertRaises(TypeError):
                syutyu.msd(path, method='custom', dt_ps=1.0)
            with self.assertRaises(ValueError) as caught:
                syutyu.msd(path, method='custom', dt_ps=1.0, coordinates='both')
            self.assertEqual(str(caught.exception), "coordinates must be 'unwrapped' or 'wrapped'")

    def test_invalid_selection_becomes_value_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = tiny_dump(Path(tmp))
            with self.assertRaises(ValueError) as caught:
                syutyu.rdf(path, method='custom', pairs=[('foo bar', 'all')], rmax=3)
            self.assertTrue(str(caught.exception).startswith('Invalid selection'))

    def test_non_orthogonal_universe_rejected(self):
        u = memory_universe([0, 10], [[[0., 0, 0], [1., 0, 0]], [[1., 0, 0], [2., 0, 0]]],
                            [[10., 10, 10]] * 2, ['1', '1'], angles=60.0)
        with self.assertRaises(ValueError) as caught:
            syutyu.rdf(u, method='custom', pairs=[('all', 'all')], rmax=3)
        self.assertEqual(str(caught.exception), 'Only orthogonal cells are supported')

    def test_positions_are_absolute(self):
        bounds = ((5, 15), (5, 15), (5, 15))
        dumps = [('id type x y z', '1 1 6.5 6 6'),
                 ('id type xs ys zs', '1 1 .15 .1 .1'),
                 ('id type xu yu zu', '1 1 6.5 6 6')]
        with tempfile.TemporaryDirectory() as tmp:
            for number, (columns, row) in enumerate(dumps):
                with self.subTest(columns=columns):
                    t = read_trajectory(write_dump(Path(tmp) / f'dump-{number}', columns,
                                                   [[row], [row]], boxes=[bounds, bounds]))
                    np.testing.assert_allclose(t.atoms.positions[0], [6.5, 6, 6], atol=1e-5)

    def test_missing_mol_and_type_columns(self):
        rows = [['1 0 0 0', '2 1 0 0']] * 2
        with tempfile.TemporaryDirectory() as tmp:
            path = write_dump(Path(tmp) / 'dump', 'id x y z', rows, steps=[0, 10])
            t = read_trajectory(path)
            np.testing.assert_array_equal(t.atoms.resids, [0, 0])
            self.assertEqual(list(t.atoms.types), ['', ''])
            self.assertEqual(list(t.atoms.elements), ['', ''])
            result = syutyu.rdf(path, method='custom', pairs=[('all', 'all')], rmax=3, bins=3,
                                intermolecular=True)
            self.assertEqual(result.tables['rdf']['CN_0'].to_numpy()[-1], 1)
            with self.assertRaises(ValueError) as caught:
                syutyu.rdf(path, method='custom', pairs=[('type 1', 'all')], rmax=3)
            self.assertEqual(str(caught.exception), 'Selection matches no atoms: type 1')
            with self.assertRaises(ValueError) as caught:
                syutyu.msd(path, method='custom', dt_ps=1.0, coordinates='wrapped',
                           com_masses={'1': 1.0})
            self.assertEqual(str(caught.exception), 'com_masses needs a type column in the dump')


if __name__ == '__main__':
    unittest.main()
