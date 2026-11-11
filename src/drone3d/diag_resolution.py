"""Resolution diagnosis for the RQ1-RQ3 study (Round 07 next step; EXPLORATORY).

Question: is the run-to-run difference in trip structure (tau_Y) caused by plateaus (many near-equal-energy routes with
different structure) or by non-convergence? Per (instance, lambda, condition): K=16 ABC runs at budget x10 and x30,
each followed by a local-search polish (single-coordinate Gaussian moves from the run's best key, 1000*D calls).
All runs and routes are stored.   Usage: python -I diag_resolution.py run | report
"""
import itertools
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from campaign import CACHE, RES, _env, instance_specs, BUDGET_MULT  # noqa: E402
from study import build_tables, structure_distance  # noqa: E402

OUT = os.path.join(CACHE, "resolution")
CONDS = ("F", "2D", "NL")
SEEDS = tuple(range(8001, 8017))     # K = 16
CASES = [(3201, 1.0), (3201, 0.0), (3205, 1.0), (3205, 0.0)]


def polish(prob, x, limit, seed):
    rng = np.random.default_rng(np.random.SeedSequence([seed, 77]))
    best_x, best_f = x.copy(), prob.evaluate(x)
    for _ in range(limit):
        y = best_x.copy()
        d = int(rng.integers(prob.D))
        y[d] = np.clip(y[d] + 0.1 * rng.standard_normal(), 0.0, 1.0)
        f = prob.evaluate(y)
        if f <= best_f:
            best_x, best_f = y, f
    return best_x


def run_case(seed_inst, lam, cond):
    from energy import Hardware
    from fast import FastProblem
    from heuristics import abc
    from instances import DEPOT_A, make_instance
    from legs import LegCost
    path = os.path.join(OUT, f"M_15_{seed_inst}_lam{lam}_{cond}.json")
    if os.path.exists(path):
        return path, "cached"
    spec = [s for s in instance_specs(["M"]) if s["n"] == 15 and s["seed"] == seed_inst][0]
    city, pluto, ras = _env()
    hw = Hardware()
    t0 = time.time()
    inst = make_instance(city, pluto, DEPOT_A, spec["n"], spec["seed"], mix=spec["mix"], lam=lam)
    leg = LegCost(city, inst.nodes, model="F", raster=ras)
    tabs = build_tables(inst, leg, hw)
    rec = dict(seed=seed_inst, lam=lam, conds={})
    for c in (cond,):
        runs = []
        for mult in (10, 30):
            for rs in SEEDS:
                prob = FastProblem(tabs[c], inst.requests, hw.q_max, hw.e_use_wh, relabel_seed=rs)
                res = abc(prob, mult * BUDGET_MULT * prob.D, seed=rs, N=20)
                e, ok, pieces = prob.feasible_energy(res.x)
                x2 = polish(prob, res.x, 1000 * prob.D, rs)
                e2, ok2, pieces2 = prob.feasible_energy(x2)
                runs.append(dict(mult=mult, seed=rs, e=e if ok else None, pieces=pieces if ok else None,
                                 e_pol=e2 if ok2 else None, pieces_pol=pieces2 if ok2 else None))
        rec["conds"][c] = runs
    rec["events"] = prob.D
    rec["wall"] = round(time.time() - t0, 1)
    os.makedirs(OUT, exist_ok=True)
    json.dump(rec, open(path + ".tmp", "w"))
    os.replace(path + ".tmp", path)
    return path, f"{rec['wall']}s"


def main_run():
    from joblib import Parallel, delayed
    t0 = time.time()
    jobs = [(s, l, c) for s, l in CASES for c in CONDS]
    gen = Parallel(n_jobs=11, return_as="generator_unordered")(delayed(run_case)(s, l, c) for s, l, c in jobs)
    for k, (p, msg) in enumerate(gen, 1):
        print(f"[res {k}/{len(jobs)}] {os.path.basename(p)} {msg} total {time.time() - t0:.0f}s", flush=True)


def report(out):
    import glob
    from instances import DEPOT_A  # noqa: F401
    L = ["# Çözünürlük tanısı (KEŞİFSEL; diag_resolution.py)", "",
         "Her (örnek, λ, koşul): K=16 ABC koşusu ×10 ve ×30 bütçe; her koşunun en iyi anahtarından yerel arama cilası (1000·D çağrı).", ""]
    rows = []
    for f in sorted(glob.glob(OUT + "/*.json")):
        r = json.load(open(f))
        n_req = 15
        # events->request map is rebuilt from the stored pieces length only for structure distance via event indices
        for c, runs in r["conds"].items():
            for key, label in (("e", "ham"), ("e_pol", "cilalı")):
                for mult in (10, 30):
                    es = np.array([x[key] for x in runs if x["mult"] == mult and x[key] is not None])
                    rows.append((r["seed"], r["lam"], c, mult, label, es))
    L += ["| örnek | λ | koşul | bütçe | sürüm | en iyi (Wh) | medyan/en iyi−1 (%) | en iyiye ≤%0,5 payı | ≤%1 payı | göreli IQR (%) |", "|---|---|---|---|---|---|---|---|---|---|"]
    pooled = {}
    for seed, lam, c, mult, label, es in rows:
        pooled.setdefault((seed, lam, c), []).append(es.min())
    for seed, lam, c, mult, label, es in rows:
        ref = min(pooled[(seed, lam, c)])
        L.append(f"| {seed} | {lam} | {c} | ×{mult} | {label} | {es.min():.1f} | {100*(np.median(es)/ref-1):.2f} | {np.mean(es <= ref*1.005):.2f} | {np.mean(es <= ref*1.01):.2f} | {100*(np.percentile(es,75)-np.percentile(es,25))/np.median(es):.2f} |")
    open(out, "w", encoding="utf8").write(chr(10).join(L) + chr(10))


def report_structure(out_md):
    """Plateau test: among runs within 1 % of the pooled best energy, how different are the trip structures?"""
    import glob
    from exact import trip_events
    from instances import DEPOT_A, make_instance
    city, pluto, _ = _env()
    L = ["", "## Plato testi: toplanmış en iyiye ≤ %1 olan koşular arasında sefer yapısı farkı (d_struct, çift medyanı)", "",
         "| örnek | λ | koşul | koşu sayısı (≤%1) | farklı yapı sayısı | d_struct medyan | enerji yayılımı (maks/min−1, %) |", "|---|---|---|---|---|---|---|"]
    for f in sorted(glob.glob(OUT + "/*.json")):
        r = json.load(open(f))
        reqs = make_instance(city, pluto, DEPOT_A, 15, r["seed"], mix=(0.4, 0.3, 0.3)).requests
        ev = trip_events(reqs, tuple(range(len(reqs))))
        for c, runs in r["conds"].items():
            cand = [(x["e_pol"], x["pieces_pol"]) for x in runs if x["e_pol"] is not None]
            best = min(e for e, _ in cand)
            near = [(e, p) for e, p in cand if e <= best * 1.01]
            ds = [structure_distance(a[1], b[1], ev, len(reqs)) for a, b in itertools.combinations(near, 2)]
            keys = {tuple(sorted(tuple(sorted(ev[k][0] for k in pc)) for pc in p)) for _, p in near}
            spread = 100 * (max(e for e, _ in near) / min(e for e, _ in near) - 1)
            L.append(f"| {r['seed']} | {r['lam']} | {c} | {len(near)} | {len(keys)} | {np.median(ds) if ds else float('nan'):.3f} | {spread:.2f} |")
    open(out_md, "a", encoding="utf8").write(chr(10).join(L) + chr(10))


if __name__ == "__main__":
    if sys.argv[1] == "run":
        main_run()
    else:
        out = os.path.join(RES, "RESOLUTION_DIAG_auto.md")
        report(out)
        report_structure(out)
        print(open(out, encoding="utf8").read())
