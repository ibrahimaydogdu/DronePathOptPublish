"""City geometry for D1 (FORMULATION_D1.md sections 2, 3, 5).

Local ENU projection (no pyproj), building obstacles from NYC footprints, service nodes
(roof / facade-level lambda scenario), obstacle-aware leg options (direct, over-fly).
Units: metres. Altitudes are absolute (NAVD88 m) = ground elevation + height above ground.
"""
from __future__ import annotations
import json
import math
from dataclasses import dataclass, field

import numpy as np
from shapely import STRtree
from shapely.geometry import LineString, MultiPolygon, Point, Polygon, shape

FT = 0.3048
R_EARTH = 6371008.8
NON_BUILDING_CODES = {"1000", "1001", "1002", "1003", "1004", "1005", "2110"}  # NYC_FIELD_DEFINITIONS.md
CEILING_AGL = 121.92  # 400 ft, 14 CFR 107.51(b)


class Projection:
    """Equirectangular projection around (lat0, lon0); error < 0.1 % over 5 km."""

    def __init__(self, lat0: float, lon0: float):
        self.lat0 = lat0
        self.lon0 = lon0
        self.cos0 = math.cos(math.radians(lat0))

    def fwd(self, lon, lat):
        x = R_EARTH * self.cos0 * np.radians(np.asarray(lon) - self.lon0)
        y = R_EARTH * np.radians(np.asarray(lat) - self.lat0)
        return x, y


def haversine(lat1, lon1, lat2, lon2):
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dphi = p2 - p1
    dl = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * R_EARTH * math.asin(math.sqrt(a))


@dataclass
class Node:
    x: float
    y: float
    z: float
    host: int = -1          # index into City.geoms, -1 = none
    kind: str = "service"   # service | depot


@dataclass
class City:
    geoms: list
    top: np.ndarray          # absolute top (m)
    ground: np.ndarray       # absolute ground (m), NaN -> filled by nearest neighbour
    height: np.ndarray       # above ground (m); 0 or NaN = unknown
    obstacle: np.ndarray     # bool, height > 0
    bin_id: list
    bbl: list
    code: list
    proj: Projection
    tree: STRtree = field(init=False, repr=False)
    obs_idx: np.ndarray = field(init=False, repr=False)

    def __post_init__(self):
        self.obs_idx = np.nonzero(self.obstacle)[0]
        self.tree = STRtree([self.geoms[i] for i in self.obs_idx])

    def terrain(self, x: float, y: float) -> float:
        if not hasattr(self, "_kd"):
            from scipy.spatial import cKDTree
            c = np.array([[g.centroid.x, g.centroid.y] for g in self.geoms])
            self._kd = cKDTree(c)
        _, i = self._kd.query([x, y])
        return float(self.ground[i])


def load_city(path: str, lat0: float = 40.725, lon0: float = -73.995) -> City:
    rows = json.load(open(path, encoding="utf8"))
    proj = Projection(lat0, lon0)
    geoms, top, ground, height, bins, bbls, codes = [], [], [], [], [], [], []
    for r in rows:
        g = shape(r["the_geom"])
        polys = []
        for p in (g.geoms if isinstance(g, MultiPolygon) else [g]):
            xs, ys = proj.fwd([c[0] for c in p.exterior.coords], [c[1] for c in p.exterior.coords])
            polys.append(Polygon(list(zip(xs, ys))).buffer(0))
        geoms.append(MultiPolygon([q for q in polys if not q.is_empty and q.geom_type == "Polygon"]) if len(polys) > 1 else polys[0])
        h = r.get("height_roof")
        gr = r.get("ground_elevation")
        height.append(float(h) * FT if h not in (None, "") else float("nan"))
        ground.append(float(gr) * FT if gr not in (None, "") else float("nan"))
        bins.append(r.get("bin"))
        bbls.append(r.get("base_bbl"))
        codes.append(r.get("feature_code"))
    height = np.array(height)
    ground = np.array(ground)
    if np.isnan(ground).any():  # fill from nearest neighbour with known ground
        cen = np.array([[g.centroid.x, g.centroid.y] for g in geoms])
        from scipy.spatial import cKDTree
        ok = ~np.isnan(ground)
        kd = cKDTree(cen[ok])
        _, j = kd.query(cen[~ok])
        ground[~ok] = ground[ok][j]
    obstacle = np.nan_to_num(height, nan=0.0) > 0
    top = ground + np.nan_to_num(height, nan=0.0)
    return City(geoms, top, ground, height, obstacle, bins, bbls, codes, proj)


FACADE_FALLBACKS: list = []   # buildings where no free ring point exists (reported by tests/campaign)


def make_node(city: City, b: int, lam: float = 1.0, c_s: float = 2.0, w: float = 5.0,
              kind: str = "service") -> Node:
    """Service node of building b. lam = 1: roof (z = top + c_s); lam < 1: facade level
    z = ground + lam*H + c_s, positioned w metres outside the footprint."""
    g = city.geoms[b]
    h = float(np.nan_to_num(city.height[b], nan=0.0))
    z = float(city.ground[b] + lam * h + c_s)
    if lam >= 1.0:
        p = g.representative_point()
        return Node(p.x, p.y, z, b, kind)
    c = g.centroid
    poly = g if g.geom_type == "Polygon" else max(g.geoms, key=lambda q: q.area)
    s = poly.exterior.project(c)
    bp = poly.exterior.interpolate(s)
    # nearest boundary point to the centroid, pushed outward along (boundary - centroid)
    d = np.array([bp.x - c.x, bp.y - c.y])
    n = np.linalg.norm(d)
    d = d / n if n > 1e-9 else np.array([1.0, 0.0])
    x0, y0 = bp.x + w * d[0], bp.y + w * d[1]
    # Round 08 fix: the pushed point must not fall inside (or hug) a neighbouring footprint.
    # Candidates: points on the w-offset ring of the host footprint, nearest to (x0, y0) first.
    ring = poly.buffer(w).exterior
    cands = [ring.interpolate(t * ring.length / 96.0) for t in range(96)]
    cands.sort(key=lambda q: (q.x - x0) ** 2 + (q.y - y0) ** 2)
    for margin in (c_s, 0.0):
        for q in cands:
            hits = city.tree.query(q.buffer(margin) if margin > 0 else q, predicate="intersects")
            if all(int(city.obs_idx[i]) == b for i in hits):
                return Node(q.x, q.y, z, b, kind)
    FACADE_FALLBACKS.append(b)
    return Node(x0, y0, z, b, kind)


def _pieces(line: LineString, geom):
    """Along-line parameter ranges (s0, s1) where line is inside geom."""
    inter = line.intersection(geom)
    if inter.is_empty:
        return []
    parts = list(inter.geoms) if hasattr(inter, "geoms") else [inter]
    out = []
    for q in parts:
        if q.geom_type == "Point":
            s = line.project(q)
            out.append((s, s))
        elif q.geom_type == "LineString":
            s_a, s_b = line.project(Point(q.coords[0])), line.project(Point(q.coords[-1]))
            out.append((min(s_a, s_b), max(s_a, s_b)))
    return out


def leg_options(city: City, a: Node, b: Node, delta: float = 10.0, w: float = 5.0):
    """Feasible segment lists for the leg a->b. Each option is a list of (ell, dz) segments.

    direct: one straight segment, clear of every non-host obstacle by delta (host: no clearance).
    over:   vertical climb to h = max(za, zb, tops within w of the line) + delta, level flight, vertical descent.
    """
    L = math.hypot(b.x - a.x, b.y - a.y)
    opts = []
    if L < 1e-9:
        dz = b.z - a.z
        return [("direct", [(0.0, dz)])]
    line = LineString([(a.x, a.y), (b.x, b.y)])
    cand = city.tree.query(line.buffer(w), predicate="intersects")
    cand_all = [int(city.obs_idx[i]) for i in cand]
    # direct
    ok = True
    for bi in cand_all:
        pcs = _pieces(line, city.geoms[bi])
        if not pcs:
            continue
        need = city.top[bi] + (0.0 if bi in (a.host, b.host) else delta)
        for s0, s1 in pcs:
            alt0 = a.z + (b.z - a.z) * s0 / L
            alt1 = a.z + (b.z - a.z) * s1 / L
            if min(alt0, alt1) < need - 1e-9:
                ok = False
                break
        if not ok:
            break
    if ok:
        opts.append(("direct", [(L, b.z - a.z)]))
    # over
    tops = [city.top[bi] for bi in cand_all]
    h = max([a.z, b.z] + [t + delta for t in tops])
    near_top = max(tops + [a.z, b.z])
    if h <= near_top + CEILING_AGL:  # 107.51(b): 400 ft above the highest structure within reach
        opts.append(("over", [(0.0, h - a.z), (L, 0.0), (0.0, b.z - h)]))
    return opts
