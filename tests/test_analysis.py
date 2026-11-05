"""Deterministic checks with analytical answers and an independent SciPy distance path."""
import tempfile
import unittest
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from syutyu import read_trajectory
from syutyu.kernels import minimum_image, pair_rdf, unwrap, msd_multi_origin
from support import frames_of, write_dump


class AnalysisTests(unittest.TestCase):
    def test_boundary_and_pair_normalization(self):
        box = np.full(3, 10.)
        xyz = np.array([[[.2, 0, 0], [9.8, 0, 0]]])
        edges = np.array([0., .5, 1.])
        g, cn, count = pair_rdf(xyz, box, edges, normalization='population')
        np.testing.assert_array_equal(count, [1, 0])
        np.testing.assert_allclose(cn, [1, 1])
        self.assertAlmostEqual(g[0], 1 / ((2 / 1000) * (4 * np.pi / 3 * .5**3)))
        np.testing.assert_allclose(minimum_image(np.array([9.6, 0, 0]), box), [-.4, 0, 0])

    def test_independent_periodic_tree(self):
        rng = np.random.default_rng(72)
        box = np.array([8., 10., 12.])
        xyz = rng.uniform(size=(1, 120, 3)) * box
        edges = np.linspace(0, 4, 21)
        _, _, count = pair_rdf(xyz, box, edges, normalization='population')
        tree = cKDTree(xyz[0], boxsize=box)
        # Tree counts ordered pairs including 120 zero-distance self-pairs.
        cumulative = tree.count_neighbors(tree, edges)
        np.testing.assert_array_equal(count, np.diff(cumulative) // 2)

    def test_ideal_gas(self):
        rng = np.random.default_rng(2)
        xyz = rng.uniform(0, 20, size=(12, 300, 3))
        edges = np.linspace(0, 10, 31)
        g, _, _ = pair_rdf(xyz, np.full(3, 20.), edges, normalization='population')
        # Large, populated shells avoid small-radius shot noise.
        self.assertLess(np.max(np.abs(g[10:] - (1 - 1/300))), .05)

    def test_radius_limit(self):
        with self.assertRaises(ValueError):
            pair_rdf(np.zeros((1, 2, 3)), np.full(3, 10), np.array([0, 6]),
                     normalization='population')

    def test_unwrap_and_ballistic_msd(self):
        t = np.arange(20.)
        real = np.zeros((20, 2, 3))
        real[:, 0, 0] = 9 + .8 * t
        real[:, 1, 0] = 1 - .8 * t
        restored = unwrap(real % 10, np.full(3, 10.))
        np.testing.assert_allclose(restored, real, atol=1e-12)
        msd, samples = msd_multi_origin(restored)
        np.testing.assert_allclose(msd, (.8 * t)**2, atol=1e-12)
        np.testing.assert_array_equal(samples, 2 * (20 - t))

    def test_pair_rdf_normalization_is_explicit(self):
        positions = np.array([[[.2, 0, 0], [9.8, 0, 0]]])
        box, edges = np.full(3, 10.), np.array([0., .5, 1.])
        with self.assertRaises(TypeError):
            pair_rdf(positions, box, edges)
        with self.assertRaises(ValueError):
            pair_rdf(positions, box, edges, normalization='bogus')
        finite = pair_rdf(positions, box, edges, normalization='finite')[0]
        population = pair_rdf(positions, box, edges, normalization='population')[0]
        np.testing.assert_allclose(finite, 2 * population)

    def test_dump_sorting(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = write_dump(Path(tmp) / 'test.dump', 'id element x y z',
                              [['2 Ar 9 0 0', '1 Ar 1 0 0'], ['1 Ar 2 0 0', '2 Ar 8 0 0']],
                              steps=[0, 100])
            trajectory = read_trajectory(path)
            steps, xyz, boxes, _velocities = frames_of(trajectory)
            np.testing.assert_array_equal(steps, [0, 100])
            np.testing.assert_array_equal(xyz[:, :, 0], [[1, 9], [2, 8]])
            np.testing.assert_array_equal(boxes[0], [10, 10, 10])


if __name__ == '__main__':
    unittest.main(verbosity=2)
