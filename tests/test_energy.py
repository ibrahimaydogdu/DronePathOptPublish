import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from energy import Hardware, segment_energy_wh, MODELS

HW = Hardware()


class TestCalibration(unittest.TestCase):
    def test_hover_endpoints(self):
        self.assertAlmostEqual(HW.hover_power(0.0), 1011.4875, places=3)
        self.assertAlmostEqual(HW.hover_power(HW.q_max), 2022.975, places=3)

    def test_manual_hover_time_reproduced(self):
        self.assertAlmostEqual(HW.e_use_wh / HW.hover_power(0.0) * 60, 32.0, places=6)
        self.assertAlmostEqual(HW.e_use_wh / HW.hover_power(HW.q_max) * 60, 16.0, places=6)

    def test_alpha_k(self):
        self.assertAlmostEqual(HW.alpha, 1.41589, places=4)
        self.assertAlmostEqual(HW.k, 41.7455, places=3)

    def test_momentum_theory_consistency(self):
        for load in (0.0, HW.q_max):
            eta = HW.ideal_hover_power(load) / HW.hover_power(load)
            self.assertTrue(0.45 < eta < 0.55, eta)

    def test_capacity(self):
        self.assertAlmostEqual(HW.q_max, 6.0)
        self.assertAlmostEqual(HW.e_use_wh, 539.46, places=6)


class TestReductions(unittest.TestCase):
    def test_level_flight_closed_form(self):
        for model in MODELS:
            if model == "S0":
                continue
            e = segment_energy_wh(model, 1000.0, 0.0, 0.0, HW, v=10.0)
            self.assertAlmostEqual(e, HW.hover_power(0.0) * 100.0 / 3600.0, places=9, msg=model)

    def test_2d_ignores_altitude(self):
        a = segment_energy_wh("2D", 500.0, 80.0, 3.0, HW)
        b = segment_energy_wh("2D", 500.0, 0.0, 3.0, HW)
        self.assertAlmostEqual(a, b, places=12)

    def test_zero_length_zero_dz(self):
        for model in MODELS:
            self.assertEqual(segment_energy_wh(model, 0.0, 0.0, 2.0, HW), 0.0)

    def test_climb_costs_more_than_level_for_F_Mo_L_S0(self):
        for model in ("F", "Mo", "L", "S0"):
            level = segment_energy_wh(model, 300.0, 0.0, 2.0, HW)
            climb = segment_energy_wh(model, 300.0, 60.0, 2.0, HW)
            self.assertGreater(climb, level, model)

    def test_D_climb_equals_level_when_horizontal_time_dominates(self):
        a = segment_energy_wh("D", 300.0, 60.0, 2.0, HW)   # 30 s horizontal, 12 s vertical
        b = segment_energy_wh("D", 300.0, 0.0, 2.0, HW)
        self.assertAlmostEqual(a, b, places=12)


class TestMonotonic(unittest.TestCase):
    def test_load_monotone(self):
        for model in MODELS:
            prev = -1.0
            for load in (0.0, 1.0, 2.5, 4.0, 6.0):
                e = segment_energy_wh(model, 400.0, 30.0, load, HW)
                self.assertGreater(e, prev, (model, load))
                prev = e

    def test_length_monotone(self):
        for model in MODELS:
            self.assertLess(segment_energy_wh(model, 100.0, 0.0, 1.0, HW),
                            segment_energy_wh(model, 200.0, 0.0, 1.0, HW))

    def test_load_bounds(self):
        with self.assertRaises(ValueError):
            segment_energy_wh("F", 100.0, 0.0, 6.5, HW)


class TestUnits(unittest.TestCase):
    def test_wh_scale(self):
        # 1 km level flight, 6 kg load: ~2023 W for 100 s = 56.2 Wh
        e = segment_energy_wh("D", 1000.0, 0.0, 6.0, HW)
        self.assertAlmostEqual(e, 56.193, places=2)

    def test_full_load_1km_within_battery(self):
        self.assertLess(segment_energy_wh("F", 1000.0, 0.0, 6.0, HW), HW.e_use_wh)


if __name__ == "__main__":
    unittest.main()
