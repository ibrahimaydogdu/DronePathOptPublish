import itertools
import os
import random
import sys
import time
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

from exact import INF, solve_exact
from fast import FastProblem, table_from_fn
from heuristics import Problem
from test_exact import random_instance

QMAX, EUSE = 4.0, 60.0


def make(seed, n_req):
    leg, reqs = random_instance(random.Random(seed), n_req)
    E = table_from_fn(leg.energy, len(leg.pts), QMAX)
    return leg, reqs, E


class TestFastDecoder(unittest.TestCase):
    def test_matches_python_decoder_when_feasible(self):
        for seed in range(1, 8):
            leg, reqs, E = make(seed, 4)
            slow = Problem(leg, reqs, QMAX, EUSE, max_trip_events=24)
            fast = FastProblem(E, reqs, QMAX, EUSE, kmax=24)
            rng = np.random.default_rng(seed)
            for _ in range(40):
                x = rng.random(slow.D)
                a = slow.evaluate(x)
                e, ok, _ = fast.feasible_energy(x)
                if a < 1e5:
                    self.assertTrue(ok)
                    self.assertAlmostEqual(a, e, places=7)
                else:
                    self.assertFalse(ok)

    def test_enumeration_reaches_exact_optimum(self):
        checked = 0
        for seed in range(1, 12):
            leg, reqs, E = make(seed, 3)
            opt, _ = solve_exact(leg, reqs, QMAX, EUSE)
            fast = FastProblem(E, reqs, QMAX, EUSE, kmax=24)
            if opt == INF or fast.D > 6:
                continue
            best = INF
            for perm in itertools.permutations(range(fast.D)):
                x = np.empty(fast.D)
                for rank, k in enumerate(perm):
                    x[k] = (rank + 0.5) / fast.D
                e, ok, _ = fast.feasible_energy(x)
                if ok and e < best:
                    best = e
            self.assertAlmostEqual(best, opt, places=7)
            checked += 1
        self.assertGreater(checked, 3)

    def test_never_below_exact(self):
        for seed in range(1, 8):
            leg, reqs, E = make(seed, 4)
            opt, _ = solve_exact(leg, reqs, QMAX, EUSE)
            fast = FastProblem(E, reqs, QMAX, EUSE)
            rng = np.random.default_rng(seed)
            for _ in range(80):
                e, ok, _ = fast.feasible_energy(rng.random(fast.D))
                if ok and opt < INF:
                    self.assertGreaterEqual(e, opt - 1e-7)

    def test_penalty_gives_gradient_and_is_larger_when_infeasible(self):
        leg, reqs, E = make(5, 4)
        fast = FastProblem(E, reqs, 4.0, 1.0)         # absurdly small battery: everything infeasible
        rng = np.random.default_rng(3)
        vals = {round(fast.evaluate(rng.random(fast.D)), 3) for _ in range(30)}
        self.assertGreater(len(vals), 5)               # not a flat penalty
        e, ok, _ = fast.feasible_energy(rng.random(fast.D))
        self.assertFalse(ok)

    def test_feasible_rate_reasonable_for_larger_n(self):
        rng0 = random.Random(7)
        leg, reqs = random_instance(rng0, 25)
        E = table_from_fn(leg.energy, len(leg.pts), QMAX)
        fast = FastProblem(E, reqs, QMAX, 400.0)
        rng = np.random.default_rng(1)
        feas = sum(fast.feasible_energy(rng.random(fast.D))[1] for _ in range(100))
        self.assertGreater(feas, 0)

    def test_speed(self):
        rng0 = random.Random(8)
        leg, reqs = random_instance(rng0, 50)
        E = table_from_fn(leg.energy, len(leg.pts), QMAX)
        fast = FastProblem(E, reqs, QMAX, 400.0)
        x = np.random.default_rng(2).random(fast.D)
        fast.evaluate(x)                                # compile
        t = time.time()
        for _ in range(200):
            fast.evaluate(x)
        per = (time.time() - t) / 200
        self.assertLess(per, 0.005)


if __name__ == "__main__":
    unittest.main()
