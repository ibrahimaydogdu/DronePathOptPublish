import math
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
from shapely.geometry import Polygon

from energy import Hardware, segment_energy_wh
from geometry import City, Node, Projection, haversine, leg_options, make_node
from legs import LegCost

HW = Hardware()


def toy_city(wall_h=50.0, ground=0.0):
    wall = Polygon([(90, -100), (110, -100), (110, 100), (90, 100)])
    host_a = Polygon([(-20, -10), (0, -10), (0, 10), (-20, 10)])
    host_b = Polygon([(200, -10), (220, -10), (220, 10), (200, 10)])
    geoms = [wall, host_a, host_b]
    height = np.array([wall_h, 20.0, 20.0])
    gr = np.full(3, ground)
    return City(geoms, gr + height, gr, height, np.ones(3, bool), [1, 2, 3], [1, 2, 3], ["2100"] * 3,
                Projection(40.725, -73.995))


class TestProjection(unittest.TestCase):
    def test_against_haversine(self):
        p = Projection(40.725, -73.995)
        for lat, lon in [(40.70, -74.02), (40.75, -73.97), (40.70, -73.97), (40.75, -74.02)]:
            x, y = p.fwd(lon, lat)
            d_proj = math.hypot(float(x), float(y))
            d_true = haversine(40.725, -73.995, lat, lon)
            self.assertLess(abs(d_proj - d_true) / d_true, 1e-3)


class TestWall(unittest.TestCase):
    def setUp(self):
        self.city = toy_city()

    def nodes(self, za=22.0, zb=22.0):
        a = Node(-10.0, 0.0, za, 1)   # on roof of host A (top 20) + 2
        b = Node(210.0, 0.0, zb, 2)
        return a, b

    def test_direct_blocked_by_wall(self):
        a, b = self.nodes()
        kinds = [k for k, _ in leg_options(self.city, a, b, delta=10.0, w=5.0)]
        self.assertNotIn("direct", kinds)
        self.assertIn("over", kinds)

    def test_over_altitude_is_wall_top_plus_delta(self):
        a, b = self.nodes()
        opts = dict(leg_options(self.city, a, b, delta=10.0, w=5.0))
        segs = opts["over"]
        self.assertAlmostEqual(segs[0][1], 60.0 - 22.0)   # climb to 50 + 10
        self.assertAlmostEqual(segs[1][0], 220.0)
        self.assertAlmostEqual(segs[2][1], 22.0 - 60.0)

    def test_direct_allowed_when_flying_high(self):
        a, b = self.nodes(za=70.0, zb=70.0)
        kinds = [k for k, _ in leg_options(self.city, a, b, delta=10.0, w=5.0)]
        self.assertIn("direct", kinds)

    def test_wall_detour_energy_matches_hand_calculation(self):
        a, b = self.nodes()
        lc = LegCost(self.city, [a, b], model="F", hw=HW)
        e = lc.energy(0, 1, 0.0)
        m = HW.m0
        ph = HW.hover_power(0.0)
        up_t = 38.0 / 5.0
        cruise_t = 22.0 * 10.0 / 10.0 * 1.0   # 220 m / 10 m/s
        down_t = 38.0 / 3.0
        e_hand = (ph * (up_t + 22.0 + down_t) + m * 9.81 * 38.0 / HW.eta_h) / 3600.0
        self.assertAlmostEqual(e, e_hand, places=9)

    def test_no_wall_equals_closed_form(self):
        city = toy_city(wall_h=0.0)
        city.obstacle[0] = False
        city.__post_init__()
        a, b = self.nodes()
        a = Node(-10.0, 0.0, 22.0, 1)
        b = Node(210.0, 0.0, 22.0, 2)
        lc = LegCost(city, [a, b], model="D", hw=HW)
        self.assertAlmostEqual(lc.energy(0, 1, 0.0), HW.hover_power(0.0) * 22.0 / 3600.0, places=9)

    def test_2d_model_same_for_all_levels(self):
        a, b = self.nodes()
        e2d = LegCost(self.city, [a, b], model="2D").energy(0, 1, 2.0)
        e2d_b = LegCost(self.city, [Node(-10.0, 0.0, 22.0, 1), Node(210.0, 0.0, 22.0, 2)], model="2D").energy(0, 1, 2.0)
        self.assertAlmostEqual(e2d, e2d_b, places=12)

    def test_load_monotone_on_leg(self):
        a, b = self.nodes()
        lc = LegCost(self.city, [a, b], model="F")
        self.assertLess(lc.energy(0, 1, 0.0), lc.energy(0, 1, 3.0))


class TestFacade(unittest.TestCase):
    def test_facade_node_is_outside_footprint(self):
        city = toy_city()
        n = make_node(city, 1, lam=0.5, c_s=2.0, w=5.0)
        self.assertFalse(city.geoms[1].contains(__import__("shapely").geometry.Point(n.x, n.y)))
        self.assertAlmostEqual(n.z, 0.0 + 0.5 * 20.0 + 2.0)

    def test_facade_node_not_inside_neighbour(self):
        from shapely.geometry import Point
        city = toy_city()
        for b in range(len(city.geoms)):
            n = make_node(city, b, lam=0.0)
            for j in range(len(city.geoms)):
                if j != b:
                    self.assertFalse(city.geoms[j].contains(Point(n.x, n.y)))

    def test_roof_node_inside_footprint(self):
        city = toy_city()
        n = make_node(city, 1, lam=1.0)
        self.assertTrue(city.geoms[1].contains(__import__("shapely").geometry.Point(n.x, n.y)))
        self.assertAlmostEqual(n.z, 22.0)


if __name__ == "__main__":
    unittest.main()
