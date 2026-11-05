"""VACF unit tests: Colab multiple-origin definition, units, parsing and rejections."""
from dataclasses import replace
from pathlib import Path
import tempfile
import unittest

import numpy as np

from syutyu import Trajectory, read_trajectory
from syutyu.analysis.vacf import velocity_autocorrelation
from support import write_dump


def velocity_dump(path, velocities, steps=None, types=None):
    frames = []
    for values in velocities:
        rows = []
        for index, velocity in enumerate(values):
            atom_type = 1 if types is None else types[index]
            rows.append(f'{index + 1} {atom_type} 0 0 0 '
                        f'{float(velocity[0])!r} {float(velocity[1])!r} {float(velocity[2])!r}')
        frames.append(rows)
    return read_trajectory(write_dump(path, 'id type x y z vx vy vz', frames, steps))


def explicit_vacf(velocities, origin_mode):
    """Independent per-atom, per-origin loop oracle returning component averages."""
    velocities = np.asarray(velocities, dtype=float)
    frames, atoms, _ = velocities.shape
    expected = np.zeros((frames, 3))
    for lag in range(frames):
        total = np.zeros(3)
        count = 0
        if origin_mode == 'multiple':
            for origin in range(frames - lag):
                for atom in range(atoms):
                    total += velocities[origin, atom] * velocities[origin + lag, atom]
                    count += 1
        else:
            for atom in range(atoms):
                total += velocities[0, atom] * velocities[lag, atom]
                count += 1
        expected[lag] = total / count
    return expected


class VelocityAutocorrelationTests(unittest.TestCase):
    def test_hand_calculated_one_atom_two_frames(self):
        with tempfile.TemporaryDirectory() as tmp:
            t = velocity_dump(Path(tmp) / 'dump', [[[1, 0, 0]], [[-1, 0, 0]]], steps=[0, 10])
            result = velocity_autocorrelation(t, ['all'], dt_ps=.001, velocity_unit='A/ps')
            table = result.tables['vacf']
            np.testing.assert_allclose(table['lag_ps'].to_numpy(), [0, .01])
            np.testing.assert_allclose(table['VACF_0_x_A2_ps2'].to_numpy(), [1, -1])
            np.testing.assert_allclose(table['VACF_0_y_A2_ps2'].to_numpy(), [0, 0])
            np.testing.assert_allclose(table['VACF_0_z_A2_ps2'].to_numpy(), [0, 0])
            np.testing.assert_allclose(table['VACF_0_A2_ps2'].to_numpy(), [1, -1])
            np.testing.assert_allclose(table['VACF_0_normalized'].to_numpy(), [1, -1])
            np.testing.assert_array_equal(table['origin_atom_pairs_0'].to_numpy(), [2, 1])
            self.assertEqual(table.columns, ['lag_ps', 'VACF_0_x_A2_ps2', 'VACF_0_y_A2_ps2',
                                             'VACF_0_z_A2_ps2', 'VACF_0_A2_ps2',
                                             'VACF_0_normalized', 'origin_atom_pairs_0'])
            self.assertEqual(result.summary['kind'], 'vacf')
            self.assertEqual(result.summary['method'], 'custom')
            self.assertEqual(result.summary['origin_mode'], 'multiple')
            self.assertEqual(result.summary['velocity_unit_input'], 'A/ps')
            self.assertEqual(result.summary['velocity_unit_output'], 'A/ps')
            self.assertAlmostEqual(result.summary['dt_ps'], .001)
            self.assertAlmostEqual(result.summary['frame_interval_ps'], .01)
            self.assertEqual(result.summary['selected_steps'], [0, 10])
            self.assertEqual(result.summary['selected_frames'], 2)
            self.assertEqual(result.summary['frame_stride'], 1)
            self.assertFalse(result.summary['drift_removed'])
            self.assertEqual(result.summary['vacf_selections'],
                             [{'selection': 'all', 'column': 'VACF_0_A2_ps2',
                               'normalized': 'VACF_0_normalized', 'atoms': 1}])
            self.assertEqual(result.summary['source_dump'], str((Path(tmp) / 'dump').resolve()))

    def test_multiple_and_single_match_independent_oracle_and_differ(self):
        with tempfile.TemporaryDirectory() as tmp:
            t = velocity_dump(Path(tmp) / 'dump', [[[1, 0, 0]], [[1, 1, 0]], [[0, -2, 0]]])
            multiple = velocity_autocorrelation(t, ['all'], dt_ps=.001, velocity_unit='A/ps')
            single = velocity_autocorrelation(t, ['all'], dt_ps=.001, velocity_unit='A/ps',
                                              origin_mode='single')
            expected_multiple = explicit_vacf(t.velocities, 'multiple')
            expected_single = explicit_vacf(t.velocities, 'single')
            table = multiple.tables['vacf']
            for axis, name in enumerate(('x', 'y', 'z')):
                np.testing.assert_allclose(table[f'VACF_0_{name}_A2_ps2'].to_numpy(),
                                           expected_multiple[:, axis])
                np.testing.assert_allclose(single.tables['vacf'][f'VACF_0_{name}_A2_ps2'].to_numpy(),
                                           expected_single[:, axis])
            np.testing.assert_allclose(table['VACF_0_A2_ps2'].to_numpy(), expected_multiple.sum(axis=1))
            np.testing.assert_allclose(single.tables['vacf']['VACF_0_A2_ps2'].to_numpy(),
                                       expected_single.sum(axis=1))
            self.assertNotAlmostEqual(table['VACF_0_A2_ps2'].to_numpy()[1],
                                      single.tables['vacf']['VACF_0_A2_ps2'].to_numpy()[1])
            np.testing.assert_array_equal(table['origin_atom_pairs_0'].to_numpy(), [3, 2, 1])
            np.testing.assert_array_equal(single.tables['vacf']['origin_atom_pairs_0'].to_numpy(),
                                          [1, 1, 1])
            self.assertEqual(single.summary['origin_mode'], 'single')
            self.assertTrue(any('Single origin' in note for note in single.summary['limitations']))
            np.testing.assert_allclose(table['VACF_0_normalized'].to_numpy(), [1, -.5 / (7 / 3), 0])

    def test_components_sum_and_no_cross_atom_correlations(self):
        with tempfile.TemporaryDirectory() as tmp:
            t = velocity_dump(Path(tmp) / 'dump',
                              [[[1, 0, 0], [0, 1, 0]], [[0, 1, 0], [1, 0, 0]]])
            table = velocity_autocorrelation(t, ['all'], dt_ps=.001,
                                             velocity_unit='A/ps').tables['vacf']
            np.testing.assert_allclose(table['VACF_0_A2_ps2'].to_numpy(), [1, 0])
            np.testing.assert_allclose(
                table['VACF_0_x_A2_ps2'].to_numpy() + table['VACF_0_y_A2_ps2'].to_numpy()
                + table['VACF_0_z_A2_ps2'].to_numpy(), table['VACF_0_A2_ps2'].to_numpy())
            np.testing.assert_allclose(table['VACF_0_x_A2_ps2'].to_numpy(), [.5, 0])
            np.testing.assert_allclose(table['VACF_0_y_A2_ps2'].to_numpy(), [.5, 0])

            general = velocity_dump(Path(tmp) / 'general', [
                [[1, -2, .5], [0, 3, 1]],
                [[.5, 1, -1], [-1, -1, .5]],
                [[-2, 0, 2], [2, 2, -3]],
            ])
            expected = explicit_vacf(general.velocities, 'multiple')
            result = velocity_autocorrelation(general, ['all'], dt_ps=.001,
                                              velocity_unit='A/ps').tables['vacf']
            for axis, name in enumerate(('x', 'y', 'z')):
                np.testing.assert_allclose(result[f'VACF_0_{name}_A2_ps2'], expected[:, axis])
            np.testing.assert_allclose(result['VACF_0_A2_ps2'].to_numpy(), expected.sum(axis=1))
            np.testing.assert_allclose(
                result['VACF_0_A2_ps2'].to_numpy(),
                result['VACF_0_x_A2_ps2'].to_numpy() + result['VACF_0_y_A2_ps2'].to_numpy()
                + result['VACF_0_z_A2_ps2'].to_numpy())
            np.testing.assert_allclose(result['VACF_0_normalized'].to_numpy()[0], 1)
            np.testing.assert_allclose(result['VACF_0_A2_ps2'].to_numpy()[0],
                                       np.mean(np.sum(general.velocities**2, axis=-1)))

    def test_constant_velocity_keeps_normalized_one(self):
        velocities = [[[1, 0, 0], [0, 1, 0]]] * 4
        with tempfile.TemporaryDirectory() as tmp:
            t = velocity_dump(Path(tmp) / 'dump', velocities)
            table = velocity_autocorrelation(t, ['all'], dt_ps=.001,
                                             velocity_unit='A/ps').tables['vacf']
            np.testing.assert_allclose(table['VACF_0_A2_ps2'].to_numpy(), [1, 1, 1, 1])
            np.testing.assert_allclose(table['VACF_0_normalized'].to_numpy(), [1, 1, 1, 1])
            np.testing.assert_array_equal(table['origin_atom_pairs_0'].to_numpy(), [8, 6, 4, 2])

    def test_velocity_unit_conversion_factor_is_one_million(self):
        with tempfile.TemporaryDirectory() as tmp:
            t = velocity_dump(Path(tmp) / 'dump', [[[2, 0, 0]], [[-1, 0, 0]]])
            ps = velocity_autocorrelation(t, ['all'], dt_ps=.001, velocity_unit='A/ps')
            fs = velocity_autocorrelation(t, ['all'], dt_ps=.001, velocity_unit='A/fs')
            np.testing.assert_allclose(fs.tables['vacf']['VACF_0_A2_ps2'].to_numpy(),
                                       ps.tables['vacf']['VACF_0_A2_ps2'].to_numpy() * 1e6)
            np.testing.assert_allclose(fs.tables['vacf']['VACF_0_normalized'].to_numpy(),
                                       ps.tables['vacf']['VACF_0_normalized'].to_numpy())
            np.testing.assert_allclose(fs.tables['vacf']['lag_ps'].to_numpy(),
                                       ps.tables['vacf']['lag_ps'].to_numpy())
            self.assertAlmostEqual(fs.summary['mean_velocity_rms_A_ps'],
                                   1000 * ps.summary['mean_velocity_rms_A_ps'])
            self.assertEqual(fs.summary['velocity_unit_input'], 'A/fs')
            self.assertEqual(fs.summary['velocity_unit_output'], 'A/ps')

    def test_reordered_dump_velocities_follow_sorted_ids(self):
        with tempfile.TemporaryDirectory() as tmp:
            path, sorted_path = Path(tmp) / 'dump', Path(tmp) / 'sorted'
            t = read_trajectory(write_dump(path, 'id type x y z vx vy vz',
                                           [['2 2 0 0 0 4 5 6', '1 1 0 0 0 1 2 3'],
                                            ['1 1 0 0 0 -1 -2 -3', '2 2 0 0 0 0 1 2']]))
            np.testing.assert_array_equal(t.ids, [1, 2])
            np.testing.assert_allclose(t.velocities[0], [[1, 2, 3], [4, 5, 6]])
            np.testing.assert_allclose(t.velocities[1], [[-1, -2, -3], [0, 1, 2]])
            sorted_t = velocity_dump(sorted_path, t.velocities)
            reordered = velocity_autocorrelation(t, ['all'], dt_ps=.001,
                                                 velocity_unit='A/ps').tables['vacf']
            ordered = velocity_autocorrelation(sorted_t, ['all'], dt_ps=.001,
                                               velocity_unit='A/ps').tables['vacf']
            for key in ordered.columns:
                np.testing.assert_allclose(reordered[key].to_numpy(), ordered[key].to_numpy())

    def test_parsing_rejects_missing_partial_and_nonfinite_velocities(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            old = read_trajectory(write_dump(tmp / 'old', 'id type x y z',
                                             [['1 1 0 0 0'], ['1 1 1 0 0']]))
            self.assertIsNone(old.velocities)
            self.assertTrue(old.select('type:1').all())
            with self.assertRaisesRegex(ValueError, 'velocities'):
                velocity_autocorrelation(old, ['all'], dt_ps=.001, velocity_unit='A/ps')
            with self.assertRaisesRegex(ValueError, 'vx/vy/vz'):
                read_trajectory(write_dump(tmp / 'partial', 'id type x y z vx',
                                           [['1 1 0 0 0 1'], ['1 1 1 0 0 1']]))
            with self.assertRaisesRegex(ValueError, 'vx/vy/vz'):
                read_trajectory(write_dump(tmp / 'partial2', 'id type x y z vx vy',
                                           [['1 1 0 0 0 1 0'], ['1 1 1 0 0 1 0']]))
            with self.assertRaisesRegex(ValueError, 'Non-finite velocities'):
                read_trajectory(write_dump(tmp / 'nan', 'id type x y z vx vy vz',
                                           [['1 1 0 0 0 nan 0 0'], ['1 1 1 0 0 1 0 0']]))
            with self.assertRaisesRegex(ValueError, 'Non-finite velocities'):
                read_trajectory(write_dump(tmp / 'inf', 'id type x y z vx vy vz',
                                           [['1 1 0 0 0 1 0 0'], ['1 1 1 0 0 inf 0 0']]))

    def test_nonfinite_velocity_on_manual_trajectory_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            t = velocity_dump(Path(tmp) / 'dump', [[[1, 0, 0]], [[-1, 0, 0]]])
            bad = replace(t, velocities=t.velocities.copy())
            bad.velocities[1, 0, 2] = np.nan
            with self.assertRaisesRegex(ValueError, 'Non-finite velocities'):
                velocity_autocorrelation(bad, ['all'], dt_ps=.001, velocity_unit='A/ps')

    def test_zero_cv0_rejected_for_undefined_normalization(self):
        with tempfile.TemporaryDirectory() as tmp:
            t = velocity_dump(Path(tmp) / 'dump', [[[0, 0, 0]], [[0, 0, 0]]])
            with self.assertRaisesRegex(ValueError, 'Cv\\(0\\) is zero'):
                velocity_autocorrelation(t, ['all'], dt_ps=.001, velocity_unit='A/ps')

    def test_stride_start_step_and_selection_counts(self):
        types = [1, 1, 2, 2]
        velocities = [[[float((m + 1) * (i + 1)), .5 * m, -float(i)] for i in range(4)]
                      for m in range(5)]
        with tempfile.TemporaryDirectory() as tmp:
            t = velocity_dump(Path(tmp) / 'dump', velocities, types=types)
            strided = velocity_autocorrelation(t, ['type:2'], dt_ps=.001,
                                               velocity_unit='A/ps', stride=2)
            table = strided.tables['vacf']
            np.testing.assert_allclose(table['lag_ps'].to_numpy(), [0, .02, .04])
            self.assertAlmostEqual(strided.summary['frame_interval_ps'], .02)
            self.assertEqual(strided.summary['selected_steps'], [0, 40])
            self.assertEqual(strided.summary['selected_frames'], 3)
            self.assertEqual(strided.summary['frame_stride'], 2)
            np.testing.assert_array_equal(table['origin_atom_pairs_0'].to_numpy(), [6, 4, 2])
            self.assertEqual(strided.summary['vacf_selections'][0]['atoms'], 2)
            expected = explicit_vacf(t.velocities[::2][:, 2:4], 'multiple')
            np.testing.assert_allclose(table['VACF_0_A2_ps2'].to_numpy(), expected.sum(axis=1))

            started = velocity_autocorrelation(t, ['all'], dt_ps=.001, velocity_unit='A/ps',
                                               start_step=20)
            np.testing.assert_allclose(started.tables['vacf']['lag_ps'].to_numpy(), [0, .01, .02])
            self.assertAlmostEqual(started.summary['frame_interval_ps'], .01)
            self.assertEqual(started.summary['selected_steps'], [20, 40])
            self.assertEqual(started.summary['selected_frames'], 3)
            np.testing.assert_array_equal(started.tables['vacf']['origin_atom_pairs_0'].to_numpy(),
                                          [12, 8, 4])
            with self.assertRaisesRegex(ValueError, 'start_step must not exceed stop_step'):
                velocity_autocorrelation(t, ['all'], dt_ps=.001, velocity_unit='A/ps',
                                         start_step=45)
            with self.assertRaisesRegex(ValueError, 'stride'):
                velocity_autocorrelation(t, ['all'], dt_ps=.001, velocity_unit='A/ps', stride=0)

    def test_max_lag_is_inclusive_and_rejects_out_of_range(self):
        velocities = [[[1, 0, 0]], [[0, 2, 0]], [[-1, 1, 0]], [[2, 0, 1]], [[0, -1, 3]]]
        with tempfile.TemporaryDirectory() as tmp:
            t = velocity_dump(Path(tmp) / 'dump', velocities)
            full = velocity_autocorrelation(t, ['all'], dt_ps=.001, velocity_unit='A/ps')
            np.testing.assert_allclose(full.tables['vacf']['lag_ps'].to_numpy(), [0, .01, .02, .03, .04])
            truncated = velocity_autocorrelation(t, ['all'], dt_ps=.001, velocity_unit='A/ps',
                                                 max_lag_ps=.02)
            np.testing.assert_allclose(truncated.tables['vacf']['lag_ps'].to_numpy(), [0, .01, .02])
            np.testing.assert_allclose(truncated.tables['vacf']['VACF_0_A2_ps2'].to_numpy(),
                                       full.tables['vacf']['VACF_0_A2_ps2'].to_numpy()[:3])
            np.testing.assert_array_equal(truncated.tables['vacf']['origin_atom_pairs_0'].to_numpy(),
                                          [5, 4, 3])
            below = velocity_autocorrelation(t, ['all'], dt_ps=.001, velocity_unit='A/ps',
                                             max_lag_ps=.019)
            np.testing.assert_allclose(below.tables['vacf']['lag_ps'].to_numpy(), [0, .01])
            np.testing.assert_array_equal(below.tables['vacf']['origin_atom_pairs_0'].to_numpy(), [5, 4])
            exact = velocity_autocorrelation(t, ['all'], dt_ps=.001, velocity_unit='A/ps',
                                             max_lag_ps=.04)
            self.assertEqual(len(exact.tables['vacf']['lag_ps'].to_numpy()), 5)
            with self.assertRaisesRegex(ValueError, 'available lag range'):
                velocity_autocorrelation(t, ['all'], dt_ps=.001, velocity_unit='A/ps',
                                         max_lag_ps=.041)
            for invalid in (0, -1, np.nan, np.inf):
                with self.assertRaisesRegex(ValueError, 'positive and finite'):
                    velocity_autocorrelation(t, ['all'], dt_ps=.001, velocity_unit='A/ps',
                                             max_lag_ps=invalid)

    def test_nonuniform_or_too_few_frames_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            nonuniform = velocity_dump(tmp / 'nonuniform',
                                       [[[1, 0, 0]], [[0, 1, 0]], [[0, 0, 1]]],
                                       steps=[0, 10, 30])
            with self.assertRaisesRegex(ValueError, 'uniformly spaced'):
                velocity_autocorrelation(nonuniform, ['all'], dt_ps=.001,
                                         velocity_unit='A/ps')
            single = velocity_dump(tmp / 'single', [[[1, 0, 0]]])
            with self.assertRaisesRegex(ValueError, 'uniformly spaced'):
                velocity_autocorrelation(single, ['all'], dt_ps=.001, velocity_unit='A/ps')
            uniform = velocity_dump(tmp / 'uniform', [[[1, 0, 0]], [[-1, 0, 0]]])
            with self.assertRaisesRegex(ValueError, 'uniformly spaced'):
                velocity_autocorrelation(uniform, ['all'], dt_ps=.001, velocity_unit='A/ps',
                                         stride=3)

    def test_argument_validation(self):
        with tempfile.TemporaryDirectory() as tmp:
            t = velocity_dump(Path(tmp) / 'dump', [[[1, 0, 0]], [[-1, 0, 0]]])
            with self.assertRaisesRegex(ValueError, 'dt_ps'):
                velocity_autocorrelation(t, ['all'], dt_ps=0, velocity_unit='A/ps')
            with self.assertRaisesRegex(ValueError, 'dt_ps'):
                velocity_autocorrelation(t, ['all'], dt_ps=np.nan, velocity_unit='A/ps')
            with self.assertRaisesRegex(ValueError, 'velocity_unit'):
                velocity_autocorrelation(t, ['all'], dt_ps=.001, velocity_unit='m/s')
            with self.assertRaisesRegex(ValueError, 'velocity_unit'):
                velocity_autocorrelation(t, ['all'], dt_ps=.001, velocity_unit='A/ps ')
            with self.assertRaisesRegex(ValueError, 'origin_mode'):
                velocity_autocorrelation(t, ['all'], dt_ps=.001, velocity_unit='A/ps',
                                         origin_mode='all')
            with self.assertRaisesRegex(ValueError, 'at least one atom selection'):
                velocity_autocorrelation(t, [], dt_ps=.001, velocity_unit='A/ps')
            with self.assertRaisesRegex(ValueError, 'no atoms'):
                velocity_autocorrelation(t, ['type:9'], dt_ps=.001, velocity_unit='A/ps')

    def test_drift_diagnostic_is_recorded_and_not_subtracted(self):
        drift = [[[.5, 0, 0], [.5, 0, 0]]] * 3
        with tempfile.TemporaryDirectory() as tmp:
            t = velocity_dump(Path(tmp) / 'dump', drift)
            result = velocity_autocorrelation(t, ['all'], dt_ps=.001, velocity_unit='A/ps')
            table = result.tables['vacf']
            self.assertAlmostEqual(result.summary['mean_velocity_rms_A_ps'], .5)
            np.testing.assert_allclose(table['VACF_0_A2_ps2'].to_numpy(), [.25, .25, .25])
            np.testing.assert_allclose(table['VACF_0_normalized'].to_numpy(), [1, 1, 1])
            self.assertFalse(result.summary['drift_removed'])
            self.assertTrue(any('drift' in note for note in result.summary['limitations']))
            self.assertTrue(any('COM' in note for note in result.summary['limitations']))

    def test_all_output_columns_finite_and_aligned(self):
        velocities = [[[1, -2, .5], [0, 3, 1], [-1, 0, 2]],
                      [[.5, 1, -1], [-1, -1, .5], [2, 2, -3]],
                      [[-2, 0, 2], [2, 2, -3], [1, -1, 1]],
                      [[0, 1, 1], [.5, -2, 1], [-1, 1, 0]]]
        with tempfile.TemporaryDirectory() as tmp:
            t = velocity_dump(Path(tmp) / 'dump', velocities, types=[1, 2, 2])
            result = velocity_autocorrelation(t, ['type:1', 'type:2'], dt_ps=.0005,
                                              velocity_unit='A/fs')
            table = result.tables['vacf']
            lag_count = table.height
            for name in table.columns:
                values = table[name].to_numpy()
                self.assertEqual(len(values), lag_count, name)
                self.assertTrue(np.isfinite(values).all(), name)
            self.assertEqual([item['selection'] for item in result.summary['vacf_selections']],
                             ['type:1', 'type:2'])
            self.assertEqual([item['column'] for item in result.summary['vacf_selections']],
                             ['VACF_0_A2_ps2', 'VACF_1_A2_ps2'])
            self.assertEqual([item['atoms'] for item in result.summary['vacf_selections']], [1, 2])
            self.assertEqual(result.tables['vacf'].columns[7:13],
                             ['VACF_1_x_A2_ps2', 'VACF_1_y_A2_ps2', 'VACF_1_z_A2_ps2',
                              'VACF_1_A2_ps2', 'VACF_1_normalized', 'origin_atom_pairs_1'])
            self.assertTrue(np.isfinite(result.summary['mean_velocity_rms_A_ps']))

    def test_trajectory_positional_construction_keeps_velocities_optional(self):
        trajectory = Trajectory(np.array([0, 10]), np.array([1]),
                                np.zeros((2, 1, 3)), np.array([[10., 10, 10]] * 2),
                                np.array([[0., 0, 0]] * 2), {'type': np.array(['1'])},
                                None, None, 'manual.dump')
        self.assertIsNone(trajectory.velocities)
        self.assertTrue(trajectory.select('type:1').all())
        self.assertTrue(trajectory.fixed_box)


if __name__ == '__main__':
    unittest.main()
