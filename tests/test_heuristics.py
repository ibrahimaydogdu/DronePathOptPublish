import os
import random
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

from exact import INF, solve_exact
from heuristics import ALGORITHMS, Problem, init_population
from test_exact import random_instance

QMAX, EUSE = 4.0, 60.0


def make(seed, n_req=4):
    rng = random.Random(seed)
    leg, reqs = random_instance(rng, n_req)
    return leg, reqs


class TestDecoder(unittest.TestCase):
    def test_decode_never_beats_exact(self):
        for seed in range(1, 8):
            leg, reqs = make(seed)
            opt, _ = solve_exact(leg, reqs, QMAX, EUSE)
            if opt == INF:
                continue
            prob = Problem(leg, reqs, QMAX, EUSE)
            rng = np.random.default_rng(seed)
            for _ in range(60):
                v = prob.evaluate(rng.random(prob.D))
                if v < 1e5:
                    self.assertGreaterEqual(v, opt - 1e-9)

    def test_precedence_repair(self):
        leg, reqs = make(3, 4)
        prob = Problem(leg, reqs, QMAX, EUSE)
        rng = np.random.default_rng(5)
        for _ in range(50):
            seq = prob.sequence(rng.random(prob.D))
            for r, kp in prob.pickup_of.items():
                self.assertLess(seq.index(kp), seq.index(kp + 1))
            self.assertEqual(sorted(seq), list(range(prob.D)))

    def test_pieces_cover_all_events_once(self):
        leg, reqs = make(4, 4)
        prob = Problem(leg, reqs, QMAX, EUSE)
        v, pieces = prob.decode(np.random.default_rng(1).random(prob.D))
        if v < 1e5:
            flat = [k for p in pieces for k in p]
            self.assertEqual(sorted(flat), list(range(prob.D)))


class TestAlgorithms(unittest.TestCase):
    def test_budget_respected_and_deterministic(self):
        leg, reqs = make(2)
        for name, fn in ALGORITHMS.items():
            p1 = Problem(leg, reqs, QMAX, EUSE)
            p2 = Problem(leg, reqs, QMAX, EUSE)
            a = fn(p1, 500, seed=1001)
            b = fn(p2, 500, seed=1001)
            self.assertLessEqual(a.evals, 500, name)
            self.assertGreaterEqual(a.evals, 495, name)
            self.assertEqual(a.best, b.best, name)

    def test_same_initial_population_per_seed(self):
        a = init_population(8, 30, 1001)
        b = init_population(8, 30, 1001)
        self.assertTrue(np.array_equal(a, b))
        self.assertFalse(np.array_equal(a, init_population(8, 30, 1002)))

    def test_reaches_exact_optimum_on_small_instances(self):
        hits = {k: 0 for k in ALGORITHMS}
        tot = 0
        for seed in range(1, 7):
            leg, reqs = make(seed, 3)
            opt, _ = solve_exact(leg, reqs, QMAX, EUSE)
            if opt == INF:
                continue
            tot += 1
            for name, fn in ALGORITHMS.items():
                res = fn(Problem(leg, reqs, QMAX, EUSE), 3000, seed=1001)
                self.assertGreaterEqual(res.best, opt - 1e-9, name)
                if res.best <= opt * (1 + 1e-6):
                    hits[name] += 1
        self.assertGreater(tot, 3)
        for name, h in hits.items():
            self.assertGreaterEqual(h / tot, 0.5, (name, hits, tot))


if __name__ == "__main__":
    unittest.main()
