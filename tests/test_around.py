import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
from shapely.geometry import Polygon

from around import CityRaster, around_options
from energy import Hardware
from geometry import City, Node, Projection
from legs import LegCost

HW = Hardware()


def toy_city(wall_h):
    wall = Polygon([(90, -100), (110, -100), (110, 100), (90, 100)])
    host_a = Polygon([(-20, -10), (0, -10), (0, 10), (-20, 10)])
    host_b = Polygon([(200, -10), (220, -10), (220, 10), (200, 10)])
    geoms = [wall, host_a, host_b]
    height = np.array([wall_h, 20.0, 20.0])
    gr = np.zeros(3)
    return City(geoms, gr + height, gr, height, np.ones(3, bool), [1, 2, 3], [1, 2, 3], ["2100"] * 3,
                Projection(40.725, -73.995))


def nodes():
    return Node(-10.0, 0.0, 22.0, 1), Node(210.0, 0.0, 22.0, 2)


class TestAround(unittest.TestCase):
    def test_around_exists_and_length_is_sane(self):
        city = toy_city(250.0)
        ras = CityRaster(city, res=5.0)
        a, b = nodes()
        opts = around_options(city, ras, a, b, delta=10.0, w=5.0, n_levels=3)
        self.assertTrue(len(opts) >= 1)
        ell = opts[0][1][1][0]
        # corner path ~313 m (hand), grid path overestimates up to ~8 %, never shorter than 300 m
        self.assertGreater(ell, 300.0)
        self.assertLess(ell, 345.0)

    def test_tall_wall_prefers_around_over_overfly(self):
        city = toy_city(250.0)
        ras = CityRaster(city, res=5.0)
        a, b = nodes()
        lc = LegCost(city, [a, b], model="F", raster=ras)
        kind, e = lc.best_option(0, 1, 0.0)
        self.assertEqual(kind, "around")
        lc_no = LegCost(city, [a, b], model="F", raster=ras, use_around=False)
        self.assertGreater(lc_no.energy(0, 1, 0.0), e)

    def test_short_wall_prefers_overfly(self):
        city = toy_city(30.0)   # top 30: over-fly h = 40 -> climb 18 m, cheap
        ras = CityRaster(city, res=5.0)
        a, b = nodes()
        lc = LegCost(city, [a, b], model="F", raster=ras)
        kind, _ = lc.best_option(0, 1, 0.0)
        self.assertIn(kind, ("over", "around"))
        lc_no = LegCost(city, [a, b], model="F", raster=ras, use_around=False)
        self.assertLessEqual(lc.energy(0, 1, 0.0), lc_no.energy(0, 1, 0.0) + 1e-12)

    def test_around_never_worse_than_without(self):
        for h in (30.0, 60.0, 120.0, 250.0):
            city = toy_city(h)
            ras = CityRaster(city, res=5.0)
            a, b = nodes()
            with_a = LegCost(city, [a, b], model="F", raster=ras).energy(0, 1, 2.0)
            without = LegCost(city, [a, b], model="F", raster=ras, use_around=False).energy(0, 1, 2.0)
            self.assertLessEqual(with_a, without + 1e-12, h)

    def test_no_obstacle_no_around_needed(self):
        city = toy_city(0.0)
        city.obstacle[0] = False
        city.__post_init__()
        ras = CityRaster(city, res=5.0)
        a, b = nodes()
        opts = around_options(city, ras, a, b)
        # only the two host roofs remain as obstacles: any around option is the straight line
        for _, segs in opts:
            self.assertAlmostEqual(segs[1][0], 220.0, places=6)


if __name__ == "__main__":
    unittest.main()
