import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np

from analysis import a12, holm, rel_iqr, struggle


def inst(meds, q, b):
    return {"medians": dict(zip(("TLBO", "GWO", "ABC"), meds)), "q": q, "b": b}


class TestStats(unittest.TestCase):
    def test_a12_extremes(self):
        self.assertEqual(a12([1, 2, 3], [4, 5, 6]), 1.0)
        self.assertEqual(a12([4, 5, 6], [1, 2, 3]), 0.0)
        self.assertEqual(a12([1, 1], [1, 1]), 0.5)

    def test_holm_monotone_and_bounded(self):
        adj = holm([0.01, 0.04, 0.03])
        self.assertAlmostEqual(adj[0], 0.03)
        self.assertTrue(all(0 <= v <= 1 for v in adj))
        self.assertEqual(sorted(adj), sorted(adj))

    def test_rel_iqr(self):
        self.assertEqual(rel_iqr([10, 10, 10, 10]), 0.0)


class TestStruggleBoolean(unittest.TestCase):
    def test_converged_tools_do_not_trigger(self):
        ml = [inst((100, 100.5, 100.2), 0.005, 0.9) for _ in range(6)]       # tools agree within 2 %
        s = [{"hit": 1.0} for _ in range(6)]
        r = struggle(ml, s)
        self.assertFalse(r["STRUGGLE"])

    def test_disagreeing_tools_with_value_over_baseline_trigger(self):
        ml = [inst((100, 112, 106), 0.01, 0.9) for _ in range(6)]            # 12 % disagreement > max(2 %, 2q)
        r = struggle(ml, [{"hit": 1.0}])
        self.assertTrue(r["STRUGGLE_ML"])
        self.assertTrue(r["STRUGGLE"])

    def test_disagreement_but_baseline_as_good_does_not_trigger(self):
        ml = [inst((100, 112, 106), 0.01, 0.4) for _ in range(6)]            # best alg not better than baseline
        r = struggle(ml, [{"hit": 1.0}])
        self.assertFalse(r["STRUGGLE_ML"])

    def test_noise_floor_protects_against_wide_run_to_run_spread(self):
        ml = [inst((100, 108, 104), 0.06, 0.9) for _ in range(6)]            # 8 % spread but q=6 % -> 2q = 12 %
        r = struggle(ml, [{"hit": 1.0}])
        self.assertFalse(r["STRUGGLE_ML"])

    def test_s_family_trigger(self):
        s = [{"hit": 0.5}, {"hit": 0.6}, {"hit": 1.0}, {"hit": 1.0}]         # 50 % of instances below 0.8 > 25 %
        r = struggle([], s)
        self.assertTrue(r["STRUGGLE_S"])
        self.assertTrue(r["STRUGGLE"])

    def test_s_family_no_trigger(self):
        s = [{"hit": 0.9}] * 8
        self.assertFalse(struggle([], s)["STRUGGLE"])

    def test_exact_boundary_values(self):
        # r exactly equal to threshold -> not strictly greater -> not unreliable
        ml = [inst((100, 102, 101), 0.0, 0.9) for _ in range(4)]
        self.assertFalse(struggle(ml, [])["STRUGGLE_ML"])


class TestStruggleInputs(unittest.TestCase):
    def _runs(self):
        import pandas as pd
        rows = []
        rng = np.random.default_rng(0)
        # M instance: ABC best (100), TLBO 112, GWO 106; RS 140, HC 120, NN 150 (single)
        for alg, mu in (("TLBO", 112), ("GWO", 106), ("ABC", 100), ("RS", 140), ("HC", 120)):
            for rs in range(1001, 1031):
                rows.append(dict(instance="M_1", family="M", algorithm=alg, run_seed=rs, feasible=True,
                                 energy=mu + rng.normal(0, 0.5), exact=""))
        rows.append(dict(instance="M_1", family="M", algorithm="NN2OPT", run_seed=0, feasible=True, energy=150.0, exact=""))
        # S instance: exact 50; ABC hits 24/30, others fewer; one infeasible ABC run
        for alg, hit in (("TLBO", 10), ("GWO", 5), ("ABC", 24)):
            for k in range(30):
                rows.append(dict(instance="S1_1", family="S1", algorithm=alg, run_seed=1001 + k,
                                 feasible=not (alg == "ABC" and k == 29),
                                 energy=50.0 if k < hit else 55.0, exact=50.0))
        return pd.DataFrame(rows)

    def test_inputs_shapes_and_values(self):
        from analysis import struggle_inputs
        ml, s = struggle_inputs(self._runs())
        self.assertEqual(len(ml), 1)
        self.assertEqual(ml[0]["best_alg"], "ABC")
        self.assertEqual(ml[0]["best_base"], "HC")
        self.assertGreater(ml[0]["b"], 0.99)                 # ABC (100) clearly below HC (120)
        self.assertLess(ml[0]["q"], 0.02)
        self.assertEqual(len(s), 1)
        self.assertAlmostEqual(s[0]["hit"], 24 / 30)         # infeasible run is not a hit (index 29 >= 24 anyway)

    def test_end_to_end_decision(self):
        from analysis import struggle, struggle_inputs
        ml, s = struggle_inputs(self._runs())
        r = struggle(ml, s)
        self.assertTrue(r["STRUGGLE_ML"])                    # 12 % disagreement, baseline clearly worse
        self.assertFalse(r["STRUGGLE_S"])                    # hit 0.8 is not < 0.8


if __name__ == "__main__":
    unittest.main()
