import math
import os
import random
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from exact import INF, Request, brute_force, solve_exact


class MatrixLeg:
    """Test double: symmetric distances, energy = d * (1 + 0.15 * load)."""

    def __init__(self, pts):
        self.pts = pts

    def energy(self, i, j, load):
        if i == j:
            return 0.0
        d = math.dist(self.pts[i], self.pts[j])
        return d * (1.0 + 0.15 * load)


def random_instance(rng, n_req, kinds="DPB"):
    n_nodes = 1 + 2 * n_req
    pts = [(0.0, 0.0)] + [(rng.uniform(-5, 5), rng.uniform(-5, 5)) for _ in range(n_nodes - 1)]
    reqs, nxt = [], 1
    for _ in range(n_req):
        k = rng.choice(kinds)
        q = round(rng.uniform(0.5, 2.5), 1)
        if k == "D":
            reqs.append(Request("D", 0, nxt, q)); nxt += 1
        elif k == "P":
            reqs.append(Request("P", nxt, 0, q)); nxt += 1
        else:
            reqs.append(Request("B", nxt, nxt + 1, q)); nxt += 2
    return MatrixLeg(pts), reqs


class TestExact(unittest.TestCase):
    def test_single_delivery_closed_form(self):
        leg = MatrixLeg([(0, 0), (3, 4)])
        reqs = [Request("D", 0, 1, 2.0)]
        tot, trips = solve_exact(leg, reqs, 6.0, 1e9)
        self.assertAlmostEqual(tot, 5.0 * (1 + 0.3) + 5.0 * 1.0)
        self.assertEqual(len(trips), 1)

    def test_single_pickup_closed_form(self):
        leg = MatrixLeg([(0, 0), (3, 4)])
        reqs = [Request("P", 1, 0, 2.0)]
        tot, _ = solve_exact(leg, reqs, 6.0, 1e9)
        self.assertAlmostEqual(tot, 5.0 * 1.0 + 5.0 * 1.3)

    def test_capacity_forces_two_trips(self):
        leg = MatrixLeg([(0, 0), (1, 0), (0, 1)])
        reqs = [Request("D", 0, 1, 4.0), Request("D", 0, 2, 4.0)]
        tot, trips = solve_exact(leg, reqs, 6.0, 1e9)
        self.assertEqual(len(trips), 2)
        self.assertAlmostEqual(tot, 2 * (1 * 1.6 + 1 * 1.0))

    def test_battery_limit_infeasible(self):
        leg = MatrixLeg([(0, 0), (10, 0)])
        reqs = [Request("D", 0, 1, 1.0)]
        tot, trips = solve_exact(leg, reqs, 6.0, 5.0)
        self.assertEqual(tot, INF)
        self.assertEqual(trips, [])

    def test_matches_brute_force_random(self):
        rng = random.Random(1001)
        for _ in range(25):
            leg, reqs = random_instance(rng, rng.choice([2, 3]))
            ev_count = sum(2 if r.kind == "B" else 1 for r in reqs)
            if ev_count > 6:
                continue
            a, _ = solve_exact(leg, reqs, 4.0, 1e9)
            b = brute_force(leg, reqs, 4.0, 1e9)
            self.assertAlmostEqual(a, b, places=9)

    def test_matches_brute_force_with_battery(self):
        rng = random.Random(1002)
        checked = 0
        for _ in range(40):
            leg, reqs = random_instance(rng, 3)
            ev_count = sum(2 if r.kind == "B" else 1 for r in reqs)
            if ev_count > 6:
                continue
            a, _ = solve_exact(leg, reqs, 4.0, 40.0)
            b = brute_force(leg, reqs, 4.0, 40.0)
            if a == INF or b == INF:
                self.assertEqual(a, b)
            else:
                self.assertAlmostEqual(a, b, places=9)
            checked += 1
        self.assertGreater(checked, 10)

    def test_precedence_pickup_before_delivery(self):
        leg = MatrixLeg([(0, 0), (5, 0), (-5, 0)])
        reqs = [Request("B", 1, 2, 1.0)]
        tot, trips = solve_exact(leg, reqs, 6.0, 1e9)
        order = trips[0][1]
        self.assertEqual([e[1] for e in order], ["p", "d"])
        self.assertAlmostEqual(tot, 5.0 + 10.0 * 1.15 + 5.0)


if __name__ == "__main__":
    unittest.main()
