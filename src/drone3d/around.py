"""Around-the-building leg option (FORMULATION_D1.md section 5, path 3 'around').

A raster of absolute roof tops (default 5 m cells). For a cruise level h the cells with
top + delta > h are blocked (dilated by w); the horizontal shortest path (8-connected,
skimage MCP_Geometric) is flown at constant altitude h after a vertical climb, then a
vertical descent. Path length is the grid length (over-estimates the Euclidean length by
up to ~8 %, conservative). Host buildings of the end nodes are exempt (no clearance);
a 2-cell disk around each end node is forced free.
"""
from __future__ import annotations
import math
import os

import numpy as np
from scipy import ndimage
from skimage.draw import polygon as draw_polygon
from skimage.graph import MCP_Geometric


def _sh() -> float:
    # RASTER_SHIFT=1: polygon vertices are moved by half a cell so that cell i (covering [i, i+1)) is tested at its centre
    return 0.5 if os.environ.get("RASTER_SHIFT", "0") == "1" else 0.0

from geometry import City, Node


class CityRaster:
    def __init__(self, city: City, res: float = 5.0, pad: float = 150.0):
        self.res = res
        xs0, ys0, xs1, ys1 = [], [], [], []
        for g in city.geoms:
            b = g.bounds
            xs0.append(b[0]); ys0.append(b[1]); xs1.append(b[2]); ys1.append(b[3])
        self.x0 = min(xs0) - pad
        self.y0 = min(ys0) - pad
        nx = int(math.ceil((max(xs1) + pad - self.x0) / res))
        ny = int(math.ceil((max(ys1) + pad - self.y0) / res))
        self.top = np.full((ny, nx), -np.inf, dtype=np.float32)
        order = np.argsort(city.top[city.obstacle])
        idxs = city.obs_idx[order]
        for bi in idxs:
            g = city.geoms[bi]
            for poly in (g.geoms if g.geom_type == "MultiPolygon" else [g]):
                cx = (np.array(poly.exterior.coords)[:, 0] - self.x0) / res - _sh()
                cy = (np.array(poly.exterior.coords)[:, 1] - self.y0) / res - _sh()
                rr, cc = draw_polygon(cy, cx, self.top.shape)
                self.top[rr, cc] = np.maximum(self.top[rr, cc], city.top[bi])
        self.shape = self.top.shape

    def cell(self, x: float, y: float):
        return int((y - self.y0) / self.res), int((x - self.x0) / self.res)


def around_options(city: City, raster: CityRaster, a: Node, b: Node, delta: float = 10.0,
                   w: float = 5.0, n_levels: int = 3):
    """List of ('around', segments) options; segments = [(0,h-za), (ell,0), (0,zb-h)]."""
    res = raster.res
    L = math.hypot(b.x - a.x, b.y - a.y)
    if L < 2 * res:
        return []
    margin = 0.35 * L + 60.0
    xa, ya, xb, yb = a.x, a.y, b.x, b.y
    ra, ca = raster.cell(xa, ya)
    rb, cb = raster.cell(xb, yb)
    r0 = max(min(ra, rb) - int(margin / res), 0)
    r1 = min(max(ra, rb) + int(margin / res) + 1, raster.shape[0])
    c0 = max(min(ca, cb) - int(margin / res), 0)
    c1 = min(max(ca, cb) + int(margin / res) + 1, raster.shape[1])
    win = raster.top[r0:r1, c0:c1]
    finite = win[np.isfinite(win)]
    h_min = max(a.z, b.z)
    top_max = float(finite.max()) if finite.size else -np.inf
    h_over = max(h_min, top_max + delta)
    if h_over <= h_min + 1e-9:
        return []
    out = []
    dil = max(int(math.ceil(w / res)), 1)
    se = ndimage.generate_binary_structure(2, 2)
    for f in np.linspace(0.0, 1.0, n_levels + 1)[:-1]:
        h = h_min + f * (h_over - h_min)
        blocked = (win + delta) > h
        # host buildings of end nodes are exempt from clearance (need top <= h)
        for nd in (a, b):
            if nd.host >= 0 and nd.host in range(len(city.geoms)) and city.top[nd.host] <= h:
                g = city.geoms[nd.host]
                for poly in (g.geoms if g.geom_type == "MultiPolygon" else [g]):
                    cx = (np.array(poly.exterior.coords)[:, 0] - raster.x0) / res - c0 - _sh()
                    cy = (np.array(poly.exterior.coords)[:, 1] - raster.y0) / res - r0 - _sh()
                    rr, cc = draw_polygon(cy, cx, blocked.shape)
                    blocked[rr, cc] = False
        raw_blocked = blocked.copy()
        blocked = ndimage.binary_dilation(blocked, structure=se, iterations=dil)
        for (r, c) in ((ra - r0, ca - c0), (rb - r0, cb - c0)):
            sq = (slice(max(r - 2, 0), r + 3), slice(max(c - 2, 0), c + 3))
            if os.environ.get("AROUND_CLEAR", "raw") == "raw":
                # Round 11 E3 fix: near an end node only the widening margin is relaxed; cells whose roof + delta is above h stay blocked.
                blocked[sq] = raw_blocked[sq]
                blocked[min(max(r, 0), blocked.shape[0] - 1), min(max(c, 0), blocked.shape[1] - 1)] = False
            else:
                blocked[sq] = False
        cost = np.where(blocked, np.inf, 1.0)
        m = MCP_Geometric(cost, fully_connected=True)
        cum, _ = m.find_costs([(ra - r0, ca - c0)], [(rb - r0, cb - c0)])
        d = cum[rb - r0, cb - c0]
        if not np.isfinite(d):
            continue
        ell = float(d) * res
        out.append(("around", [(0.0, h - a.z), (ell, 0.0), (0.0, b.z - h)]))
    return out
