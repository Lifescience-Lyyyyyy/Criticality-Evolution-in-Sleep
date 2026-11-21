import unittest

import numpy as np

from observation import continuous_signal, rate_counts, spike_counts


class ObservationTest(unittest.TestCase):
    def test_half_open_edges_and_segment_boundaries(self):
        a = spike_counts([(0, .11, ([0, .039, .040, .08, .102], [])),
                          (1, 1.08, ([1., 1.04], [1.04]))])
        np.testing.assert_array_equal(a.population(), [[2, 1], [1, 2]])
        np.testing.assert_array_equal(a.rebin(2).population(), [[3], [3]])
        with self.assertRaises(ValueError):
            spike_counts([(0, .08, ([.08],))])
        with self.assertRaises(ValueError):
            spike_counts([(0, .08, ([0],)), (.07, .18, ([],))])

    def test_rate_exact_integral_and_replay(self):
        obs = rate_counts([np.array([10., 20., 30., 40.])], .025, units=3, seed=19)
        # Two full bins: .025*10 + .025*20 + .025*30 + .005*40 = 1.7 counts/unit.
        self.assertTrue(np.allclose(obs.metadata["expected_unit_hz"], [21.25]*3))
        ref = rate_counts([np.ones((2000, 2))*5], .04, seed=77)
        alt = rate_counts([np.ones((2000, 2))*5], .04, seed=77)
        gamma = rate_counts([np.ones((2000, 2))*5], .04, seed=77, shared_shape=1)
        np.testing.assert_array_equal(ref.segments[0], alt.segments[0])
        self.assertFalse(np.array_equal(ref.segments[0], gamma.segments[0]))
        with self.assertRaises(ValueError):
            rate_counts([np.array([-1.])], .04, units=1, seed=77)

    def test_continuous_keeps_units_and_has_no_count_api(self):
        signal = continuous_signal([[.01, -.04], [.1]], 512, "mV")
        self.assertEqual(signal.unit, "mV")
        self.assertEqual(signal.sample_hz, 512)
        self.assertFalse(hasattr(signal, "population"))


if __name__ == "__main__":
    unittest.main()
