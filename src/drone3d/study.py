"""RQ1-RQ3 study: route decisions under alternative energy-model forms, service heights and counterfactuals.

Conditions per (instance, lambda): energy models F, D, Mo, L, S0, 2D (RQ1), NL = F with load-independent leg energy
(RQ3 counterfactual), FREE = F with straight legs and no building obstacles (RQ2 counterfactual).
Per condition: K independent ABC runs; best route = lowest feasible own-model energy; cross-evaluation of every
route under every other condition's table; same-condition split-half noise floor (PARAMETRIC_PLAN section 5, Round 06 protocol).

Usage: python -I study.py run [family n_count]   |   python -I study.py smoke
"""
from __future__ import annotations
import itertools
import json
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from campaign import CACHE, RES, _env, instance_specs, BUDGET_MULT  # noqa: E402

MODELS = ("F", "D", "Mo", "L", "S0", "2D")
CONDS = MODELS + ("NL", "FREE")
LAMS = (1.0, 0.75, 0.5, 0.0)            # 0.0 = street level
RUN_SEEDS = tuple(range(7001, 7009))     # K = 8 runs per condition
OUT = os.path.join(CACHE, os.environ.get("STUDY_OUT_MAIN", "study"))


class FreeLeg:
    """Obstacle-free leg: one straight segment (horizontal distance, altitude difference)."""
    def __init__(self, nodes, hw, model="F", v=10.0, kappa_c=1.0):
        self.nodes, self.hw, self.model, self.v, self.kappa_c = nodes, hw, model, v, kappa_c

    def options(self, i, j):
        a, b = self.nodes[i], self.nodes[j]
        return [("free", [(math.hypot(b.x - a.x, b.y - a.y), b.z - a.z)])]


def eval_route(E, reqs, events, pieces, q_max, e_use):
    """Energy (Wh) and battery/capacity feasibility of a FIXED route (list of trips, each a list of event indices)."""
    from fast import STEP
    qmax_i = int(round(q_max / STEP))
    total, ok = 0.0, True
    for piece in pieces:
        load0 = sum(int(round(reqs[events[k][0]].q / STEP)) for k in piece if reqs[events[k][0]].kind == "D")
        load, pos, e = load0, 0, 0.0
        if load > qmax_i:
            ok = False
        for k in piece:
            e += E[pos, events[k][2], min(max(load, 0), qmax_i)]
            load += int(round(events[k][3] / STEP))
            if load > qmax_i or load < 0:
                ok = False
            pos = events[k][2]
        e += E[pos, 0, min(max(load, 0), qmax_i)]
        if e > e_use + 1e-9:
            ok = False
        total += e
    return float(total), bool(ok)


def structure_distance(pa, pb, events, n_req):
    """Share of request pairs whose same-trip / different-trip status differs between two routes (0 = same trip structure)."""
    def trip_of(pieces):
        t = {}
        for ti, piece in enumerate(pieces):
            for k in piece:
                t[events[k][0]] = ti
        return t
    ta, tb = trip_of(pa), trip_of(pb)
    pairs = list(itertools.combinations(range(n_req), 2))
    diff = sum((ta[i] == ta[j]) != (tb[i] == tb[j]) for i, j in pairs)
    return diff / len(pairs) if pairs else 0.0


def order_distance(pa, pb):
    """Normalised Kendall-tau distance between the concatenated event orders (events present in both)."""
    sa = [k for p in pa for k in p]
    sb = [k for p in pb for k in p]
    pos = {k: i for i, k in enumerate(sb)}
    seq = [pos[k] for k in sa]
    n = len(seq)
    inv = sum(1 for i in range(n) for j in range(i + 1, n) if seq[i] > seq[j])
    return inv / (n * (n - 1) / 2) if n > 1 else 0.0


def build_tables(inst, nodes_leg, hw, models=MODELS, mo_variants=False):
    """Energy tables per model; NL/FREE from F (default) and, if mo_variants, NLMo/FREEMo from the Mo model (Round 10 main model)."""
    from fast import build_table_sym
    n = len(inst.nodes)
    tabs = {}
    for m in models:
        nodes_leg.model = m
        tabs[m] = build_table_sym(nodes_leg, n, hw.q_max)
    if "F" in tabs:
        tabs["NL"] = np.repeat(tabs["F"][:, :, :1], tabs["F"].shape[2], axis=2)
        tabs["FREE"] = build_table_sym(FreeLeg(inst.nodes, hw), n, hw.q_max)
    if mo_variants:
        tabs["NLMo"] = np.repeat(tabs["Mo"][:, :, :1], tabs["Mo"].shape[2], axis=2)
        tabs["FREEMo"] = build_table_sym(FreeLeg(inst.nodes, hw, model="Mo"), n, hw.q_max)
    return tabs


_REGION = {}


def _env_region(region):
    """Case A: campaign environment. Case B (Morgan P&DC): footprints main + strip N merged (derived file, Evidence untouched)."""
    if region == "A":
        return _env()
    if "B" not in _REGION:
        from around import CityRaster
        from geometry import load_city
        from instances import load_pluto
        ev = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "Evidence", "data"))
        der = os.path.join(os.path.dirname(__file__), "derived")
        os.makedirs(der, exist_ok=True)
        merged = os.path.join(der, "caseB_footprints_merged.json")
        if not os.path.exists(merged):
            a = json.load(open(os.path.join(ev, "NYC_Building_Footprints_raw.json"), encoding="utf8"))
            b = json.load(open(os.path.join(ev, "NYC_Building_Footprints_strip_N_raw.json"), encoding="utf8"))
            seen, out = set(), []
            for r in a + b:
                k = r.get("bin")
                if k in seen:
                    continue
                seen.add(k); out.append(r)
            json.dump(out, open(merged, "w", encoding="utf8"))
        pl = {}
        for f in ("NYC_MapPLUTO_raw.json", "NYC_MapPLUTO_strip_N_raw.json"):
            pl.update(load_pluto(os.path.join(ev, f)))
        city = load_city(merged)
        _REGION["B"] = (city, pl, CityRaster(city, res=5.0))
    return _REGION["B"]


def _depot_b(city):
    """Case B depot = representative point of the Morgan P&DC footprint (BIN 1012831; the OSM node is an address point
    outside the footprint, DEPOT_SOURCES.md)."""
    k = [i for i, b in enumerate(city.bin_id) if str(b) == "1012831"]
    if not k:
        raise ValueError("Morgan P&DC footprint (BIN 1012831) not found")
    g = city.geoms[k[0]]
    p = g.representative_point()
    proj = city.proj
    lon = proj.lon0 + math.degrees(p.x / (6371008.8 * proj.cos0))
    lat = proj.lat0 + math.degrees(p.y / 6371008.8)
    return {"name": "Morgan P&DC (footprint BIN 1012831)", "lat": lat, "lon": lon}


def run_cell(spec, lam, budget_mult, seeds, out_dir=OUT, conds=CONDS):
    from energy import Hardware
    from fast import FastProblem
    from heuristics import abc
    from instances import DEPOT_A, make_instance
    from legs import LegCost
    iid = f"{spec['family']}_{spec['n']}_{spec['seed']}_lam{lam}"
    path = os.path.join(out_dir, iid + ".json")
    if os.path.exists(path):
        return iid, "cached"
    region = spec.get("region", "A")
    city, pluto, ras = _env_region(region)
    hw = Hardware()
    t0 = time.time()
    inst = make_instance(city, pluto, DEPOT_A if region == "A" else _depot_b(city), spec["n"], spec["seed"], mix=spec["mix"], lam=lam)
    leg = LegCost(city, inst.nodes, model="F", raster=ras)
    tabs = build_tables(inst, leg, hw, models=tuple(m for m in MODELS if m in conds), mo_variants=("NLMo" in conds))
    t_tab = time.time() - t0
    N_ABC = 20
    cells = {}
    for c in conds:
        E = tabs[c]
        runs = []
        for rs in seeds:
            prob = FastProblem(E, inst.requests, hw.q_max, hw.e_use_wh, relabel_seed=rs)
            res = abc(prob, budget_mult * prob.D, seed=rs, N=N_ABC)
            e, ok, pieces = prob.feasible_energy(res.x)
            runs.append(dict(seed=rs, energy=e if ok else None, feasible=bool(ok), pieces=pieces))
        events = prob.events
        feas = [r for r in runs if r["feasible"]]
        best = min(feas, key=lambda r: r["energy"]) if feas else None
        if best is not None:   # internal consistency: re-evaluating the stored route under its own table must match
            e2, ok2 = eval_route(E, inst.requests, events, best["pieces"], hw.q_max, hw.e_use_wh)
            assert ok2 and abs(e2 - best["energy"]) < 1e-6 * max(1.0, e2), (iid, c, e2, best["energy"])
        cells[c] = dict(runs=runs, best=best)
    # cross evaluation + noise floor
    cross, noise = {}, {}
    n_req = len(inst.requests)
    for a in conds:
        ba = cells[a]["best"]
        for b in conds:
            bb = cells[b]["best"]
            if ba is None or bb is None:
                continue
            e_ab, ok_ab = eval_route(tabs[b], inst.requests, events, ba["pieces"], hw.q_max, hw.e_use_wh)
            cross[f"{a}>{b}"] = dict(energy=e_ab, feasible=ok_ab, own_best=bb["energy"],
                                     regret=e_ab / bb["energy"] - 1.0,
                                     d_struct=structure_distance(ba["pieces"], bb["pieces"], events, n_req),
                                     d_order=order_distance(ba["pieces"], bb["pieces"]),
                                     trips_a=len(ba["pieces"]), trips_b=len(bb["pieces"]))
        feas = [r for r in cells[a]["runs"] if r["feasible"]]
        if len(feas) >= 4:
            h = len(feas) // 2
            ha = min(feas[:h], key=lambda r: r["energy"])
            hb = min(feas[h:2 * h], key=lambda r: r["energy"])
            noise[a] = dict(rel_energy=abs(ha["energy"] - hb["energy"]) / min(ha["energy"], hb["energy"]),
                            d_struct=structure_distance(ha["pieces"], hb["pieces"], events, n_req),
                            d_order=order_distance(ha["pieces"], hb["pieces"]))
    rec = dict(instance=iid, lam=lam, n=spec["n"], seed=spec["seed"], events=len(events), budget_mult=budget_mult,
               table_seconds=round(t_tab, 1), wall=round(time.time() - t0, 1),
               own={c: (cells[c]["best"]["energy"] if cells[c]["best"] else None) for c in conds},
               feasible_runs={c: sum(r["feasible"] for r in cells[c]["runs"]) for c in conds},
               run_energies={c: [r["energy"] for r in cells[c]["runs"]] for c in conds},
               routes={c: (cells[c]["best"]["pieces"] if cells[c]["best"] else None) for c in conds},
               cross=cross, noise=noise)
    os.makedirs(out_dir, exist_ok=True)
    json.dump(rec, open(path + ".tmp", "w"))
    os.replace(path + ".tmp", path)
    return iid, f"{time.time() - t0:.0f}s"


def main_run(n_inst=8):
    from joblib import Parallel, delayed
    specs = [s for s in instance_specs(["M"]) if s["n"] == 15][:n_inst]
    jobs = [(s, lam) for s in specs for lam in LAMS]
    t0 = time.time()
    _throttle()
    gen = Parallel(n_jobs=int(os.environ.get("STUDY_JOBS", "3")), return_as="generator_unordered")(
        delayed(run_cell)(s, lam, 10 * BUDGET_MULT, RUN_SEEDS) for s, lam in jobs)
    for k, (iid, msg) in enumerate(gen, 1):
        print(f"[study {k}/{len(jobs)}] {iid} {msg} total {time.time() - t0:.0f}s", flush=True)


def main_run3(n_each=10):
    """Round 10: main model = Mo. Conditions Mo, 2D, NLMo, FREEMo on the same 28 instances x 4 lambdas, K=8, x10, same seeds."""
    from joblib import Parallel, delayed
    _throttle()
    mix0 = (0.4, 0.3, 0.3)
    base = [s for s in instance_specs(["M"]) if s["n"] == 15][:8]
    specs = base + [dict(family="RA", region="A", n=15, mix=mix0, seed=3401 + k) for k in range(n_each)] +         [dict(family="RB", region="B", n=15, mix=mix0, seed=3501 + k) for k in range(n_each)]
    jobs = [(s, lam) for s in specs for lam in LAMS]
    out = os.path.join(CACHE, os.environ.get("STUDY_OUT_DIR", "study_mo"))
    t0 = time.time()
    gen = Parallel(n_jobs=int(os.environ.get("STUDY_JOBS", "3")), return_as="generator_unordered")(
        delayed(run_cell)(s, lam, 10 * BUDGET_MULT, RUN_SEEDS, out, ("Mo", "2D", "NLMo", "FREEMo")) for s, lam in jobs)
    for k, (iid, msg) in enumerate(gen, 1):
        print(f"[study3 {k}/{len(jobs)}] {iid} {msg} total {time.time() - t0:.0f}s", flush=True)


def _throttle():
    """Low priority (Windows BELOW_NORMAL; inherited by workers) so the machine stays usable."""
    try:
        import ctypes
        ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(), 0x4000)
    except Exception:
        pass


def main_run2(n_each=10):
    """Extension: 10 new case A instances (seeds 3401+) and 10 case B instances (seeds 3501+), M n=15, 4 lambdas."""
    from joblib import Parallel, delayed
    _throttle()
    mix = (0.4, 0.3, 0.3)
    specs = [dict(family="RA", region="A", n=15, mix=mix, seed=3401 + k) for k in range(n_each)] +             [dict(family="RB", region="B", n=15, mix=mix, seed=3501 + k) for k in range(n_each)]
    jobs = [(s, lam) for s in specs for lam in LAMS]
    t0 = time.time()
    gen = Parallel(n_jobs=int(os.environ.get("STUDY_JOBS", "3")), return_as="generator_unordered")(
        delayed(run_cell)(s, lam, 10 * BUDGET_MULT, RUN_SEEDS) for s, lam in jobs)
    for k, (iid, msg) in enumerate(gen, 1):
        print(f"[study2 {k}/{len(jobs)}] {iid} {msg} total {time.time() - t0:.0f}s", flush=True)


def smoke():
    spec = dict(family="S1", n=5, mix=(0.4, 0.3, 0.3), seed=3001)
    out = os.path.join(RES, "study_smoke")
    print(run_cell(spec, 0.5, 20, (7001, 7002, 7003, 7004), out_dir=out))


if __name__ == "__main__":
    if sys.argv[1] == "smoke":
        smoke()
    elif sys.argv[1] == "run3":
        main_run3()
    elif sys.argv[1] == "run2":
        main_run2()
    elif sys.argv[1] == "smokeB":
        spec = dict(family="RB", region="B", n=5, mix=(0.4, 0.3, 0.3), seed=3501)
        print(run_cell(spec, 0.5, 20, (7001, 7002, 7003, 7004), out_dir=os.path.join(RES, "study_smokeB")))
    else:
        main_run(int(sys.argv[2]) if len(sys.argv) > 2 else 8)
