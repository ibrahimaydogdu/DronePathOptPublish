"""Campaign runner (CAMPAIGN_DESIGN.md v3, ANALYSIS_PLAN.md v2).

python -I campaign.py tune          # N tuning (instances 2001-2005, runs 2101-2110)
python -I campaign.py run [fam..]   # main campaign (families S1 S2 M L), resumable
python -I campaign.py long          # 10x budget reference on M and L (ABC, TLBO; seeds 5001-5005)
python -I campaign.py one <fam> <n> <idx>   # smoke test of one instance (prints, writes nothing)
"""
from __future__ import annotations
import csv
import hashlib
import json
import os
import platform
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
RES = os.path.join(ROOT, "04_Results")
CACHE = os.path.join(RES, "cache")
DATA = os.path.join(ROOT, "Evidence", "data")

BUDGET_MULT = 500
RUN_SEEDS = list(range(1001, 1031))
MAIN = ("TLBO", "GWO", "ABC")
BASE_STOCH = ("RS", "HC")
CHECK_MULTS = (100, 250, 500)

FAMILIES = {
    "S1": dict(ns=[6, 8], mix=(0.4, 0.3, 0.3), count=10, seed0=3001, exact=True),
    "S2": dict(ns=[5, 6], mix=(0.1, 0.1, 0.8), count=10, seed0=3101, exact=True),
    "M": dict(ns=[15, 25], mix=(0.4, 0.3, 0.3), count=10, seed0=3201, exact=False),
    "L": dict(ns=[50], mix=(0.4, 0.3, 0.3), count=10, seed0=3301, exact=False),
}

_G = {}


def _env():
    if "city" not in _G:
        from around import CityRaster
        from geometry import load_city
        from instances import load_pluto
        _G["city"] = load_city(os.path.join(DATA, "NYC_Building_Footprints_raw.json"))
        _G["pluto"] = load_pluto(os.path.join(DATA, "NYC_MapPLUTO_raw.json"))
        _G["ras"] = CityRaster(_G["city"], res=5.0)
    return _G["city"], _G["pluto"], _G["ras"]


def instance_specs(families):
    out = []
    for fam in families:
        f = FAMILIES[fam]
        for ni, n in enumerate(f["ns"]):
            for k in range(f["count"]):
                out.append(dict(family=fam, n=n, mix=f["mix"], seed=f["seed0"] + 20 * ni + k, exact=f["exact"]))
    return out


def build_instance(spec):
    import numpy as np
    from energy import Hardware
    from fast import build_table_sym
    from instances import DEPOT_A, make_instance
    from legs import LegCost
    city, pluto, ras = _env()
    hw = Hardware()
    inst = make_instance(city, pluto, DEPOT_A, spec["n"], spec["seed"], mix=spec["mix"])
    leg = LegCost(city, inst.nodes, model="F", raster=ras)
    os.makedirs(CACHE, exist_ok=True)
    key = f"{spec['family']}_{spec['n']}_{spec['seed']}.npy"
    path = os.path.join(CACHE, key)
    t0 = time.time()
    stats_path = path.replace(".npy", ".json")
    if os.path.exists(path) and os.path.exists(stats_path):
        E = np.load(path)
        built = False
    else:
        E = build_table_sym(leg, len(inst.nodes), hw.q_max)
        kinds = {}
        n_nodes = len(inst.nodes)
        for i in range(n_nodes):
            for j in range(i + 1, n_nodes):
                k, _ = leg.best_option(i, j, 2.0)
                kinds[k] = kinds.get(k, 0) + 1
        np.save(path + ".tmp.npy", E)
        os.replace(path + ".tmp.npy", path)
        json.dump({"pairs": n_nodes * (n_nodes - 1) // 2, "best_option_at_2kg": kinds,
                   "table_seconds": round(time.time() - t0, 1)}, open(stats_path, "w"))
        built = True
    return inst, leg, E, hw, (time.time() - t0 if built else 0.0)


def run_instance(spec, algos=None, run_seeds=None, budget_mult=BUDGET_MULT, Ns=None, tag="main"):
    """Return (rows, feature_row) for one instance."""
    import numpy as np
    from exact import solve_exact
    from fast import FastProblem
    from heuristics import ALGORITHMS, BASELINES
    inst, leg, E, hw, t_table = build_instance(spec)
    run_seeds = run_seeds or RUN_SEEDS
    allalg = {**ALGORITHMS, **BASELINES}
    algos = algos or list(MAIN + BASE_STOCH + ("NN2OPT",))
    Ns = Ns or _load_N()
    D = FastProblem(E, inst.requests, hw.q_max, hw.e_use_wh).D
    budget = budget_mult * D
    iid = f"{spec['family']}_{spec['n']}_{spec['seed']}"
    exact_e, exact_trips, exact_t = None, None, None
    if spec["exact"]:
        t0 = time.time()
        exact_e, trips = solve_exact(leg, inst.requests, hw.q_max, hw.e_use_wh)
        exact_t, exact_trips = time.time() - t0, len(trips)
    rows, best_known, best_pieces = [], float("inf"), None
    for alg in algos:
        fn = allalg[alg]
        seeds = [0] if alg == "NN2OPT" else run_seeds
        for rs in seeds:
            prob = FastProblem(E, inst.requests, hw.q_max, hw.e_use_wh,
                               relabel_seed=(None if alg == "NN2OPT" else rs))   # NN2OPT is deterministic: no relabelling
            t0 = time.time()
            res = fn(prob, budget, seed=rs, N=Ns.get(alg, 30)) if (alg in MAIN or alg in Ns) else fn(prob, budget, seed=rs)
            wall = time.time() - t0
            e, ok, pieces = prob.feasible_energy(res.x)
            if ok and e < best_known:
                best_known, best_pieces = e, pieces
            chk = {}
            for m in CHECK_MULTS:
                lim = m * D
                vals = [v for (u, v) in res.history if u <= lim]
                chk[m] = vals[-1] if vals else float("nan")
            rows.append(dict(instance=iid, family=spec["family"], n=spec["n"], algorithm=alg, run_seed=rs,
                             feasible=bool(ok), energy=(e if ok else float("nan")), penalised=res.best, evals=res.evals,
                             wall=wall, tie_fraction=prob.tie_fraction(res.x),
                             scouts=(res.extra or {}).get("scouts", ""), chk100=chk[100], chk250=chk[250],
                             chk500=chk[500], exact=exact_e if exact_e is not None else "", budget=budget, tag=tag))
    n_trips = exact_trips if exact_trips is not None else (len(best_pieces) if best_pieces else "")
    nB = sum(1 for r in inst.requests if r.kind == "B")
    feat = dict(instance=iid, family=spec["family"], n=spec["n"], seed=spec["seed"], events=D, B_share=nB / spec["n"],
                total_load=round(sum(r.q for r in inst.requests), 2), best_known_trips=n_trips,
                nodes=len(inst.nodes), table_seconds=round(t_table, 1), exact=exact_e if exact_e is not None else "",
                exact_seconds=round(exact_t, 1) if exact_t else "", best_known=best_known,
                option_share=json.dumps(_option_share(spec)))
    return rows, feat


def _option_share(spec):
    p = os.path.join(CACHE, f"{spec['family']}_{spec['n']}_{spec['seed']}.json")
    if not os.path.exists(p):
        return {}
    d = json.load(open(p))
    return {k: round(v / d["pairs"], 3) for k, v in d["best_option_at_2kg"].items()}


def ensure_tables(specs):
    """Build every distance table once, one job per instance (no concurrent writes to the same file)."""
    from joblib import Parallel, delayed

    def one(s):
        build_instance(s)
        return s["seed"]
    Parallel(n_jobs=11)(delayed(one)(s) for s in specs)


def _load_N():
    p = os.path.join(RES, "tuning_N.json")
    return json.load(open(p)) if os.path.exists(p) else {}


def _append(path, dict_rows):
    if not dict_rows:
        return
    new = not os.path.exists(path)
    with open(path, "a", newline="", encoding="utf8") as f:
        w = csv.DictWriter(f, fieldnames=list(dict_rows[0].keys()))
        if new:
            w.writeheader()
        w.writerows(dict_rows)


def _done(path, col="instance"):
    if not os.path.exists(path):
        return set()
    with open(path, encoding="utf8") as f:
        return {r[col] for r in csv.DictReader(f)}


def write_run_log():
    import numba
    import numpy
    import scipy
    import shapely
    import skimage
    files = {}
    for fn in sorted(os.listdir(HERE)):
        if fn.endswith(".py"):
            files[fn] = hashlib.sha256(open(os.path.join(HERE, fn), "rb").read()).hexdigest()
    data = {}
    for fn in ("NYC_Building_Footprints_raw.json", "NYC_MapPLUTO_raw.json"):
        data[fn] = hashlib.sha256(open(os.path.join(DATA, fn), "rb").read()).hexdigest()
    import joblib
    log = dict(time=time.strftime("%Y-%m-%d %H:%M:%S"), python=platform.python_version(), numpy=numpy.__version__,
               scipy=scipy.__version__, shapely=shapely.__version__, skimage=skimage.__version__, numba=numba.__version__,
               code_sha256=files, data_sha256=data, budget_mult=BUDGET_MULT, run_seeds=[RUN_SEEDS[0], RUN_SEEDS[-1]],
               decoder=dict(W_BATT=20.0, W_CAP=300.0, W_PAIR=500.0, kmax=24), joblib=joblib.__version__,
               tuning_N=_load_N(), argv=sys.argv[1:])
    name = "run_log.json" if sys.argv[1] == "run" else f"run_log_{sys.argv[1]}.json"
    json.dump(log, open(os.path.join(RES, name), "w"), indent=1)


def cmd_run(families):
    from joblib import Parallel, delayed
    write_run_log()
    specs = instance_specs(families)
    out = os.path.join(RES, "campaign_runs.csv")
    feats = os.path.join(RES, "campaign_instances.csv")
    done = _done(out)
    todo = [s for s in specs if f"{s['family']}_{s['n']}_{s['seed']}" not in done]
    print(f"{len(specs)} instances, {len(todo)} to do", flush=True)
    t0 = time.time()
    ensure_tables(todo)
    print(f"tables ready {time.time() - t0:.0f}s", flush=True)
    gen = Parallel(n_jobs=11, return_as="generator_unordered")(delayed(run_instance)(s) for s in todo)
    for k, (rows, feat) in enumerate(gen, 1):
        _append(out, rows)
        _append(feats, [feat])
        print(f"[{k}/{len(todo)}] {feat['instance']} done, elapsed {time.time() - t0:.0f}s", flush=True)


def cmd_long():
    from joblib import Parallel, delayed
    write_run_log()
    specs = instance_specs(["M", "L"])
    out = os.path.join(RES, "campaign_long.csv")
    done = _done(out)
    todo = [s for s in specs if f"{s['family']}_{s['n']}_{s['seed']}" not in done]
    gen = Parallel(n_jobs=11, return_as="generator_unordered")(
        delayed(run_instance)(s, ["ABC", "TLBO"], list(range(5001, 5006)), 10 * BUDGET_MULT, None, "long") for s in todo)
    t0 = time.time()
    for k, (rows, feat) in enumerate(gen, 1):
        _append(out, rows)
        print(f"[long {k}/{len(todo)}] {feat['instance']} {time.time() - t0:.0f}s", flush=True)


def cmd_diag():
    """Round 06 diagnostics (tag diag_*): N sweeps, GWO/HC long runs, ABC x30. Resumable per (instance, config)."""
    from joblib import Parallel, delayed
    write_run_log()
    specs = instance_specs(["M", "L"])
    m25 = [s for s in specs if s["family"] == "M" and s["n"] == 25]
    l50 = [s for s in specs if s["family"] == "L"][:4]
    cfgs = []   # (spec, alg, N, mult)
    for s in m25 + l50:
        for N in (10, 15, 50, 100):
            cfgs.append((s, "TLBO", N, 10))
        for N in (20, 100):
            cfgs.append((s, "GWO", N, 10))
        cfgs.append((s, "HC", 30, 10))
    for s in m25:
        cfgs.append((s, "ABC", 20, 30))
    out = os.path.join(RES, "campaign_diag.csv")
    done = set()
    if os.path.exists(out):
        import csv
        with open(out, encoding="utf8") as f:
            done = {(r["instance"], r["tag"]) for r in csv.DictReader(f)}
    def key(c):
        s, a, N, m = c
        return f"{s['family']}_{s['n']}_{s['seed']}", f"diag_{a}_N{N}_x{m}"
    todo = [c for c in cfgs if key(c) not in done]
    def job(c):
        s, a, N, m = c
        rows, feat = run_instance(s, [a], list(range(6001, 6006)), m * BUDGET_MULT, {a: N}, key(c)[1])
        return rows, key(c)
    gen = Parallel(n_jobs=11, return_as="generator_unordered")(delayed(job)(c) for c in todo)
    t0 = time.time()
    for k, (rows, kk) in enumerate(gen, 1):
        _append(out, rows)
        print(f"[diag {k}/{len(todo)}] {kk} {time.time() - t0:.0f}s", flush=True)


def cmd_ablate():
    """Operator-swap ablation (Round 06): ABC with all-coordinate move, TLBO with one-coordinate move; 1x budget, 10 seeds."""
    from joblib import Parallel, delayed
    import functools
    import heuristics as H
    write_run_log()
    specs = [s for s in instance_specs(["M", "L"]) if (s["family"] == "M" and s["n"] == 25) or s["family"] == "L"]
    variants = {"ABCall": functools.partial(H.abc, coords="all"), "TLBO1": functools.partial(H.tlbo, coords=1)}
    base = {"ABCall": "ABC", "TLBO1": "TLBO"}
    out = os.path.join(RES, "campaign_ablate.csv")
    done = _done(out)
    todo = [s for s in specs if f"{s['family']}_{s['n']}_{s['seed']}" not in done]

    def job(spec):
        import heuristics
        for name, fn in variants.items():
            heuristics.ALGORITHMS[name] = fn
        heuristics.ALGORITHMS_EXTRA = True
        rows, feat = run_instance(spec, list(variants), list(range(6101, 6111)), BUDGET_MULT, {k: _load_N()[base[k]] for k in variants}, "ablate")
        return rows, feat
    gen = Parallel(n_jobs=11, return_as="generator_unordered")(delayed(job)(s) for s in todo)
    t0 = time.time()
    for k, (rows, feat) in enumerate(gen, 1):
        _append(out, rows)
        print(f"[ablate {k}/{len(todo)}] {feat['instance']} {time.time() - t0:.0f}s", flush=True)


def cmd_tune():
    import numpy as np
    from joblib import Parallel, delayed
    tune = [dict(family="T", n=8, mix=(0.4, 0.3, 0.3), seed=2001, exact=False),
            dict(family="T", n=6, mix=(0.1, 0.1, 0.8), seed=2002, exact=False),
            dict(family="T", n=25, mix=(0.4, 0.3, 0.3), seed=2003, exact=False),
            dict(family="T", n=25, mix=(0.4, 0.3, 0.3), seed=2004, exact=False),
            dict(family="T", n=50, mix=(0.4, 0.3, 0.3), seed=2005, exact=False)]
    cfgs = [(a, N) for a in MAIN for N in (20, 30, 50)]
    write_run_log()
    ensure_tables(tune)

    def one(spec, alg, N):
        rows, _ = run_instance(spec, [alg], list(range(2101, 2111)), BUDGET_MULT, {alg: N}, f"tune_{alg}_{N}")
        for r in rows:
            r["N"] = N
        return rows

    res = Parallel(n_jobs=11)(delayed(one)(s, a, N) for s in tune for (a, N) in cfgs)
    flat = [r for rows in res for r in rows]
    _append(os.path.join(RES, "tuning_runs.csv"), flat)
    import pandas as pd
    df = pd.DataFrame(flat)
    df = df[df["feasible"]]
    med = df.groupby(["instance", "algorithm", "N"])["energy"].median().reset_index()
    best = {}
    table = []
    for alg in MAIN:
        sub = med[med["algorithm"] == alg].pivot(index="instance", columns="N", values="energy")
        ranks = sub.rank(axis=1, method="average")
        mr = ranks.mean(axis=0)
        table.append((alg, mr.to_dict()))
        best[alg] = int(min(mr.index, key=lambda N: (mr[N], N)))     # ties -> smaller N
    json.dump(best, open(os.path.join(RES, "tuning_N.json"), "w"))
    json.dump({a: {str(k): float(v) for k, v in t.items()} for a, t in table}, open(os.path.join(RES, "tuning_ranks.json"), "w"), indent=1)
    print("tuned N:", best)
    print("mean ranks:", table)


def cmd_one(fam, n, idx):
    spec = [s for s in instance_specs([fam]) if s["n"] == int(n)][int(idx)]
    t0 = time.time()
    rows, feat = run_instance(spec, run_seeds=list(range(1001, 1004)))
    print(feat)
    import collections
    d = collections.defaultdict(list)
    for r in rows:
        d[r["algorithm"]].append((r["energy"], r["feasible"]))
    for a, v in d.items():
        print(a, [(round(e, 1) if e == e else None, f) for e, f in v])
    print("total", round(time.time() - t0, 1), "s")


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "run":
        cmd_run(sys.argv[2:] or list(FAMILIES))
    elif cmd == "long":
        cmd_long()
    elif cmd == "ablate":
        cmd_ablate()
    elif cmd == "diag":
        cmd_diag()
    elif cmd == "tune":
        cmd_tune()
    elif cmd == "one":
        cmd_one(*sys.argv[2:5])
