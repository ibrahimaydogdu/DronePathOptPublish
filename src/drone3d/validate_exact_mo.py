"""Exact-solver gate for the Mo study (AGENTS rule 1c; Round 07/10): small instances, all four study tables.

For S-size instances (n_req requests, <= max_events events) and every lambda, the exact solver (trip enumeration + set
partition DP, solve_exact on a table-backed leg) and ABC (K=8 seeds 7001-7008, budget ratio as in the study) are run on the
Mo, 2D, NLMo and FREEMo tables.  Reported: ABC best-of-8 gap to exact, and the cross-regret of every condition's route in
the Mo table computed once with EXACT routes and once with ABC routes (does ABC reproduce the exact sensitivity?).
Usage: python -I validate_exact_mo.py [n_req] [n_instances]  -> 04_Results/VALIDATION_EXACT_MO_auto.md
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from campaign import BUDGET_MULT, RES  # noqa: E402

LAMS = (1.0, 0.75, 0.5, 0.0)
CONDS = ("Mo", "2D", "NLMo", "FREEMo")
SEEDS = tuple(range(7001, 7009))


class TableLeg:
    """Table-backed leg: energy(i, j, load_kg) from E[i, j, load index]."""

    def __init__(self, E, q_max):
        from fast import STEP
        self.E, self.step, self.nq = E, STEP, E.shape[2] - 1

    def energy(self, i, j, load_kg):
        if i == j:
            return 0.0
        return float(self.E[i, j, min(max(int(round(load_kg / self.step)), 0), self.nq)])


def run_cell(n_req, seed, lam, max_events):
    from energy import Hardware
    from exact import solve_exact, trip_events
    from fast import FastProblem
    from heuristics import abc
    from instances import DEPOT_A, make_instance
    from legs import LegCost
    from study import _env_region, build_tables, eval_route
    city, pluto, ras = _env_region("A")
    hw = Hardware()
    inst = make_instance(city, pluto, DEPOT_A, n_req, seed, mix=(0.4, 0.3, 0.3), lam=lam)
    leg = LegCost(city, inst.nodes, model="F", raster=ras)
    tabs = build_tables(inst, leg, hw, models=("Mo", "2D"), mo_variants=True)
    reqs = inst.requests
    events = trip_events(reqs, tuple(range(len(reqs))))
    idx = {(e[0], e[1]): k for k, e in enumerate(events)}
    if len(events) > max_events:
        return dict(skip=True, seed=seed, lam=lam, D=len(events))
    out = dict(skip=False, seed=seed, lam=lam, D=len(events), exact={}, abc={}, abc_med={}, ex_routes={}, abc_routes={})
    for c in CONDS:
        E = tabs[c]
        opt, trips = solve_exact(TableLeg(E, hw.q_max), reqs, hw.q_max, hw.e_use_wh)
        out["exact"][c] = opt
        out["ex_routes"][c] = [[idx[(r, t)] for (r, t, _n) in order] for _cost, order in trips]
        runs, best = [], (float("inf"), None)
        for rs in SEEDS:
            prob = FastProblem(E, reqs, hw.q_max, hw.e_use_wh, relabel_seed=rs)
            res = abc(prob, 10 * BUDGET_MULT * prob.D, seed=rs, N=20)
            e, ok, pieces = prob.feasible_energy(res.x)
            if ok:
                runs.append(e)
                if e < best[0]:
                    best = (e, pieces)
        out["abc"][c] = min(runs) if runs else None
        out["abc_med"][c] = float(np.median(runs)) if runs else None
        out["abc_routes"][c] = best[1]
    Em = tabs["Mo"]
    for c in CONDS:
        out.setdefault("cross_exact", {})[c] = eval_route(Em, reqs, events, out["ex_routes"][c], hw.q_max, hw.e_use_wh)
        out.setdefault("cross_abc", {})[c] = eval_route(Em, reqs, events, out["abc_routes"][c], hw.q_max, hw.e_use_wh) if out["abc_routes"][c] else (None, False)
    return out


def main():
    from joblib import Parallel, delayed
    n_req = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    n_inst = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    max_events = 10
    jobs = [(n_req, 3001 + k, lam) for k in range(n_inst) for lam in LAMS]
    res = Parallel(n_jobs=int(os.environ.get("STUDY_JOBS", "6")))(delayed(run_cell)(n, s, l, max_events) for n, s, l in jobs)
    json.dump(res, open(os.path.join(RES, "validate_exact_mo_v2_raw.json"), "w"))
    ok = [r for r in res if not r["skip"]]
    L = ["# Kesin çözücü kapısı, Mo çalışması (otomatik; validate_exact_mo.py; KEŞİFSEL)", "",
         f"S boyutu: n_req={n_req}, karışım (0,4/0,3/0,3), tohumlar 3001..{3000 + n_inst}, 4 λ; olay sayısı ≤ {max_events}. Hücre: {len(ok)} kullanılan / {len(res)} (atlanan: {sum(r['skip'] for r in res)}). Kesin = sefer numaralandırma + küme bölme DP (aynı tablo); ABC K=8 (7001–7008), bütçe 10×. Olay sayısı aralığı: {min(r['D'] for r in ok)}–{max(r['D'] for r in ok)}.", "",
         "## ABC en iyi-8 / kesin − 1", "| koşul | hücre | ABC=kesin (≤1e-6) payı | medyan boşluk % | maks boşluk % | koşu-medyanı boşluğu medyan % |", "|---|---|---|---|---|---|"]
    for c in CONDS:
        g = [(r["abc"][c] / r["exact"][c] - 1) for r in ok if r["abc"][c] is not None and np.isfinite(r["exact"][c])]
        gm = [(r["abc_med"][c] / r["exact"][c] - 1) for r in ok if r["abc_med"][c] is not None and np.isfinite(r["exact"][c])]
        g = np.array(g)
        L.append(f"| {c} | {len(g)} | {np.mean(g <= 1e-6):.2f} | {100*np.median(g):.2f} | {100*g.max():.2f} | {100*np.median(gm):.2f} |")
    L += ["", "## Mo tablosunda çapraz pişmanlık: kesin rotalarla ve ABC rotalarıyla (Mo kesin optimumuna göre)", "| λ | çift | n (uygun) | kesin: medyan % | ABC: medyan % | |fark| medyan (puan) | işaret uyumu (kesin>0 ⇔ ABC>0) |", "|---|---|---|---|---|---|---|"]
    for lam in LAMS:
        rs = [r for r in ok if r["lam"] == lam]
        for c in ("FREEMo", "NLMo", "2D"):
            ex, ab = [], []
            for r in rs:
                (ee, eok), (ae, aok) = r["cross_exact"][c], r["cross_abc"][c]
                if eok and aok and np.isfinite(r["exact"]["Mo"]):
                    ex.append(ee / r["exact"]["Mo"] - 1); ab.append(ae / r["exact"]["Mo"] - 1)
            if not ex:
                continue
            ex, ab = np.array(ex), np.array(ab)
            L.append(f"| {lam} | {c}→Mo | {len(ex)}/{len(rs)} | {100*np.median(ex):.1f} | {100*np.median(ab):.1f} | {100*np.median(np.abs(ex-ab)):.2f} | {np.mean((ex > 1e-9) == (ab > 1e-9)):.2f} |")
    L += ["", "Not: kesin rota Mo'da uygunsuzsa hücre dışlanır (AŞAĞI yanlı). ABC'nin Mo boşluğu çalışmadaki çözücü çözünürlüğünün S boyutundaki alt sınır tahminidir; n=15 için geçerlilik KANITI DEĞİL (n=15'te kesin çözüm yok)."]
    open(os.path.join(RES, "VALIDATION_EXACT_MO_v2_auto.md"), "w", encoding="utf8").write(chr(10).join(L) + chr(10))
    print(chr(10).join(L))


if __name__ == "__main__":
    main()

