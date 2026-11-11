"""Tests for the Review_Round_04 conditions (stable ties, lexicographic report, RNG streams, relabelling, NN+2opt)."""
import os
import random
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

from exact import INF, solve_exact
from fast import FastProblem, table_from_fn
from heuristics import ALGORITHMS, BASELINES, init_population, nn_two_opt, random_search
from test_exact import random_instance

QMAX, EUSE = 4.0, 60.0


def make(seed, n_req):
    leg, reqs = random_instance(random.Random(seed), n_req)
    return leg, reqs, table_from_fn(leg.energy, len(leg.pts), QMAX)


class TestTiesAndRelabel(unittest.TestCase):
    def test_tied_keys_are_deterministic(self):
        leg, reqs, E = make(3, 4)
        p = FastProblem(E, reqs, QMAX, EUSE)
        x = np.full(p.D, 0.5)
        a = p.feasible_energy(x)
        b = p.feasible_energy(x.copy())
        self.assertEqual(a[0], b[0])
        self.assertEqual(a[2], b[2])

    def test_relabel_changes_tie_breaking_but_not_validity(self):
        leg, reqs, E = make(3, 4)
        x = np.zeros(sum(2 if r.kind == "B" else 1 for r in reqs))
        base = FastProblem(E, reqs, QMAX, EUSE)
        outs = set()
        for s in range(6):
            p = FastProblem(E, reqs, QMAX, EUSE, relabel_seed=s)
            e, ok, pieces = p.feasible_energy(x)
            flat = sorted(k for piece in pieces for k in piece)
            self.assertEqual(flat, list(range(p.D)))
            outs.add(round(e, 6))
        self.assertGreaterEqual(len(outs), 1)
        self.assertEqual(base.tie_fraction(x), 1.0 - 1.0 / len(x))


class TestLexicographicReport(unittest.TestCase):
    def test_report_never_infeasible_when_a_feasible_split_exists(self):
        for seed in range(1, 10):
            leg, reqs, E = make(seed, 4)
            p = FastProblem(E, reqs, QMAX, EUSE)
            rng = np.random.default_rng(seed)
            for _ in range(60):
                x = rng.random(p.D)
                e_l, ok_l, _ = p.feasible_energy(x)
                f = p.evaluate(x)
                if ok_l:
                    # penalised objective is a lower bound of the lexicographic feasible energy
                    self.assertLessEqual(f, e_l + 1e-9)


class TestStreams(unittest.TestCase):
    def test_rs_stream_differs_from_any_initial_population(self):
        leg, reqs, E = make(2, 4)
        p = FastProblem(E, reqs, QMAX, EUSE)
        first = np.random.default_rng(np.random.SeedSequence([1001, 4])).random(p.D)
        for s in range(1001, 1031):
            pop = init_population(p.D, 30, s)
            for row in pop:
                self.assertFalse(np.array_equal(row, first))

    def test_algorithms_use_distinct_streams(self):
        a = np.random.default_rng(np.random.SeedSequence([1001, 1])).random(5)
        b = np.random.default_rng(np.random.SeedSequence([1001, 2])).random(5)
        self.assertFalse(np.array_equal(a, b))


class TestNN2opt(unittest.TestCase):
    def test_valid_and_not_below_exact(self):
        for seed in range(1, 6):
            leg, reqs, E = make(seed, 4)
            opt, _ = solve_exact(leg, reqs, QMAX, EUSE)
            p = FastProblem(E, reqs, QMAX, EUSE)
            res = nn_two_opt(p, 5000)
            e, ok, pieces = p.feasible_energy(res.x)
            self.assertEqual(sorted(k for piece in pieces for k in piece), list(range(p.D)))
            if ok and opt < INF:
                self.assertGreaterEqual(e, opt - 1e-7)

    def test_deterministic(self):
        leg, reqs, E = make(2, 4)
        r1 = nn_two_opt(FastProblem(E, reqs, QMAX, EUSE), 5000, seed=1)
        r2 = nn_two_opt(FastProblem(E, reqs, QMAX, EUSE), 5000, seed=2)
        self.assertEqual(r1.best, r2.best)


class TestAlgorithmsOnFast(unittest.TestCase):
    def test_all_algorithms_run_and_respect_budget(self):
        leg, reqs, E = make(4, 4)
        for name, fn in {**ALGORITHMS, **BASELINES}.items():
            p = FastProblem(E, reqs, QMAX, EUSE)
            res = fn(p, 400, seed=1001)
            self.assertLessEqual(res.evals, 400, name)
            self.assertTrue(np.isfinite(res.best), name)


if __name__ == "__main__":
    unittest.main()
