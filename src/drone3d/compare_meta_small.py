"""Small three-metaheuristic check (AGENTS rule 4): ABC, TLBO, GWO vs exact, Mo table only, same 48 small cells.
Reuses the exact optima stored by validate_exact_mo.py (validate_exact_mo_v2_raw.json). Equal budget, 8 seeds 7001-7008.
Usage: python -I compare_meta_small.py [n_req] [n_instances] -> 04_Results/COMPARE_META_SMALL_auto.md"""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from campaign import BUDGET_MULT, RES  # noqa: E402
from validate_exact_mo import LAMS, SEEDS  # noqa: E402

def run_cell(n_req, seed, lam, exact_opt):
    from energy import Hardware
    from exact import trip_events
    from fast import FastProblem
    import heuristics as H
    from instances import DEPOT_A, make_instance
    from legs import LegCost
    from study import _env_region, build_tables
    city, pluto, ras = _env_region("A")
    hw = Hardware()
    inst = make_instance(city, pluto, DEPOT_A, n_req, seed, mix=(0.4, 0.3, 0.3), lam=lam)
    leg = LegCost(city, inst.nodes, model="F", raster=ras)
    tabs = build_tables(inst, leg, hw, models=("Mo",), mo_variants=True)
    E, reqs = tabs["Mo"], inst.requests
    out = dict(seed=seed, lam=lam, exact=exact_opt)
    for name in ("ABC", "TLBO", "GWO"):
        f = H.ALGORITHMS[name]
        runs = []
        for rs in SEEDS:
            prob = FastProblem(E, reqs, hw.q_max, hw.e_use_wh, relabel_seed=rs)
            res = f(prob, 10 * BUDGET_MULT * prob.D, seed=rs, N=20)
            e, ok, _ = prob.feasible_energy(res.x)
            runs.append(e if ok else None)
        out[name] = runs
    return out

def main():
    from joblib import Parallel, delayed
    n_req = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    n_inst = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    stored = json.load(open(os.path.join(RES, "validate_exact_mo_v2_raw.json")))
    ex = {(r["seed"], r["lam"]): r["exact"]["Mo"] for r in stored if not r["skip"]}
    jobs = [(3001 + k, lam) for k in range(n_inst) for lam in LAMS if (3001 + k, lam) in ex]
    res = Parallel(n_jobs=int(os.environ.get("STUDY_JOBS", "12")))(delayed(run_cell)(n_req, s, l, ex[(s, l)]) for s, l in jobs)
    json.dump(res, open(os.path.join(RES, "compare_meta_small_raw.json"), "w"))
    L = ["# Üç meta-sezgiselin küçük kıyası (otomatik; compare_meta_small.py; KEŞİFSEL)", "",
         f"Mo tablosu, {len(res)} hücre (n_req={n_req}, tohum 3001..{3000+n_inst}, 4 λ), kesin optimum doğrulama dosyasından; her yöntem 8 koşu (7001–7008), N=20, eşit değerlendirme bütçesi.", "",
         "| yöntem | hücre | en iyi-8 = kesin (≤1e-6) | medyan boşluk % (en iyi-8) | maks boşluk % | koşu-medyanı boşluğu medyan % | koşu başına kesin payı | uygunsuz koşu |", "|---|---|---|---|---|---|---|---|"]
    for a in ("ABC", "TLBO", "GWO"):
        best, med, hit, bad, tot = [], [], [], 0, 0
        for r in res:
            v = [x for x in r[a] if x is not None]
            bad += len(r[a]) - len(v); tot += len(r[a])
            if not v or not np.isfinite(r["exact"]): continue
            best.append(min(v) / r["exact"] - 1); med.append(np.median(v) / r["exact"] - 1)
            hit += [x / r["exact"] - 1 <= 1e-6 for x in v]
        b = np.array(best)
        L.append(f"| {a} | {len(b)} | {np.mean(b <= 1e-6):.2f} | {100*np.median(b):.2f} | {100*b.max():.2f} | {100*np.median(med):.2f} | {np.mean(hit):.2f} | {bad}/{tot} |")
    open(os.path.join(RES, "COMPARE_META_SMALL_auto.md"), "w", encoding="utf8").write("\n".join(L) + "\n")
    print("\n".join(L))

if __name__ == "__main__":
    main()
