"""Instance generator for D1 (FORMULATION_D1.md section 3).

Geometry is real (NYC footprints, MapPLUTO). Demand is SYNTHETIC: customer buildings are sampled
from PLUTO-joined buildings within a radius of the depot with weight 1 + unitstotal [assumption],
request kinds and masses follow user-set parameters [assumption].
"""
from __future__ import annotations
import json
import math
import random
from dataclasses import dataclass, field

import numpy as np
from shapely.geometry import Point

from exact import Request
from geometry import City, Node, make_node

DEPOT_A = {"name": "Church Street Station", "lat": 40.7127418, "lon": -74.0106511}   # DATA_SOURCES.md (OSM way/277904718)
DEPOT_B = {"name": "Morgan P&DC", "lat": 40.7506257, "lon": -73.9985278}             # DATA_SOURCES.md (OSM node/2703228136)
ELIGIBLE_LANDUSE = {"1", "2", "3", "4", "5"}   # PLUTO land use codes: residential .. commercial/office


def load_pluto(path: str) -> dict:
    return {str(int(float(r["bbl"]))): r for r in json.load(open(path, encoding="utf8"))}


def host_building(city: City, x: float, y: float) -> int:
    p = Point(x, y)
    hits = city.tree.query(p, predicate="intersects")
    cands = [int(city.obs_idx[i]) for i in hits if city.geoms[int(city.obs_idx[i])].contains(p)]
    if not cands:
        raise ValueError("depot point is not inside any footprint of the loaded data")
    return max(cands, key=lambda b: city.top[b])


@dataclass
class Instance:
    nodes: list[Node]
    requests: list[Request]
    meta: dict = field(default_factory=dict)


def eligible_buildings(city: City, pluto: dict, x0: float, y0: float, radius: float):
    out, weights = [], []
    for i in range(len(city.geoms)):
        if not city.obstacle[i] or city.code[i] != "2100":
            continue
        c = city.geoms[i].centroid
        if math.hypot(c.x - x0, c.y - y0) > radius:
            continue
        try:
            key = str(int(float(city.bbl[i])))
        except (TypeError, ValueError):
            continue
        row = pluto.get(key)
        if row is None or row.get("landuse") not in ELIGIBLE_LANDUSE:
            continue
        units = float(row.get("unitstotal") or 0.0)
        out.append(i)
        weights.append(1.0 + units)
    return out, np.array(weights)


def make_instance(city: City, pluto: dict, depot: dict, n_req: int, seed: int, radius: float = 1500.0,
                  mix: tuple = (0.4, 0.3, 0.3), q_range: tuple = (0.5, 3.0), lam: float = 1.0) -> Instance:
    rng = random.Random(seed)
    xs, ys = city.proj.fwd(depot["lon"], depot["lat"])
    x0, y0 = float(xs), float(ys)
    hb = host_building(city, x0, y0)
    depot_node = make_node(city, hb, lam=1.0, kind="depot")
    cand_all, w_all = eligible_buildings(city, pluto, x0, y0, radius)
    keep = [k for k, c in enumerate(cand_all) if c != hb]
    cand = [cand_all[k] for k in keep]
    w = [float(w_all[k]) for k in keep]
    nodes = [depot_node]
    reqs: list[Request] = []

    def draw_building(exclude=()):
        for _ in range(1000):
            b = rng.choices(cand, weights=w, k=1)[0]
            if b not in exclude:
                return b
        raise RuntimeError("could not draw a building")

    kinds = rng.choices(["D", "P", "B"], weights=mix, k=n_req)
    for kind in kinds:
        q = round(rng.uniform(*q_range), 1)
        if kind == "B":
            b1 = draw_building()
            b2 = draw_building(exclude=(b1,))
            nodes.append(make_node(city, b1, lam=lam))
            nodes.append(make_node(city, b2, lam=lam))
            reqs.append(Request("B", len(nodes) - 2, len(nodes) - 1, q))
        else:
            b1 = draw_building()
            nodes.append(make_node(city, b1, lam=lam))
            if kind == "D":
                reqs.append(Request("D", 0, len(nodes) - 1, q))
            else:
                reqs.append(Request("P", len(nodes) - 1, 0, q))
    meta = {"depot": depot["name"], "depot_building_bin": city.bin_id[hb], "radius_m": radius, "seed": seed,
            "n_req": n_req, "mix": mix, "q_range": q_range, "lam": lam, "demand": "SYNTHETIC (PLUTO-weighted)"}
    return Instance(nodes, reqs, meta)
