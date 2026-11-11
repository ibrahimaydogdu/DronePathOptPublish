"""Horizontal / vertical decomposition of the obstacle effect (Round 09 open item; EXPLORATORY).

For the stored best F route of every cell (cache/study/*.json) each leg (at its actual load) is re-evaluated under the
F model with its best geometric option, and compared with the obstacle-free straight leg at the SAME node altitudes
(the FREE condition).  Energy terms (model F):
    H = P_h * ell / v                      level-flight energy over the horizontal path length
    T = P_h * (t - ell / v)                extra hover-power time spent climbing/descending (when vertical time dominates)
    C = m g dz+ / eta_h                    climb (potential-energy) term
Excess = actual - straight; additive per leg: E_F(route) - E_FREE(route) = dH + dT + dC.
Checks: sum of actual legs == stored own F energy, and sum of straight legs == stored cross F>FREE energy.
Usage: python -I decompose.py   -> 04_Results/DECOMP_auto.md
"""
import glob
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from campaign import CACHE, RES, instance_specs  # noqa: E402
from energy import G, Hardware  # noqa: E402
from exact import trip_events  # noqa: E402
from fast import STEP  # noqa: E402

V = 10.0


def terms(ell, dz, load_kg, hw):
    up, down = max(dz, 0.0), max(-dz, 0.0)
    m = hw.m0 + load_kg
    ph = hw.hover_power(load_kg)
    t = max(ell / V, up / hw.v_up, down / hw.v_down)
    return ph * (ell / V) / 3600.0, ph * (t - ell / V) / 3600.0, m * G * up / hw.eta_h / 3600.0


def leg_terms(segs, load_kg, hw):
    H = T = C = 0.0
    for ell, dz in segs:
        h, t, c = terms(ell, dz, load_kg, hw)
        H, T, C = H + h, T + t, C + c
    return np.array([H, T, C])


def decompose_cell(path):
    from instances import DEPOT_A, make_instance
    from legs import LegCost
    from study import _depot_b, _env_region
    rec = json.load(open(path))
    fam, n, seed = rec["instance"].split("_lam")[0].split("_")[0], rec["n"], rec["seed"]
    lam = rec["lam"]
    region = "B" if fam == "RB" else "A"
    mix = (0.4, 0.3, 0.3)
    if fam == "M":
        mix = [s for s in instance_specs(["M"]) if s["n"] == 15 and s["seed"] == seed][0]["mix"]
    city, pluto, ras = _env_region(region)
    hw = Hardware()
    inst = make_instance(city, pluto, DEPOT_A if region == "A" else _depot_b(city), n, seed, mix=mix, lam=lam)
    leg = LegCost(city, inst.nodes, model="F", raster=ras)
    reqs = inst.requests
    events = trip_events(reqs, tuple(range(len(reqs))))
    out = dict(act=np.zeros(3), free=np.zeros(3), kinds={"direct": 0, "over": 0, "around": 0, "free": 0, "other": 0},
               path_act=0.0, path_free=0.0, up_act=0.0, up_free=0.0, e_act=0.0, e_free=0.0)

    def leg_add(i, j, load_kg):
        opts = leg.options(i, j)
        if i == j:
            return
        best = None
        for kind, segs in opts:
            e = sum(sum(terms(ell, dz, load_kg, hw)) for ell, dz in segs)
            if best is None or e < best[0]:
                best = (e, kind, segs)
        e, kind, segs = best
        a, b = inst.nodes[i], inst.nodes[j]
        L = math.hypot(b.x - a.x, b.y - a.y)
        free = leg_terms([(L, b.z - a.z)], load_kg, hw)
        act = leg_terms(segs, load_kg, hw)
        out["act"] += act
        out["free"] += free
        out["kinds"][kind if kind in out["kinds"] else "other"] += 1
        out["path_act"] += sum(s[0] for s in segs)
        out["path_free"] += L
        out["up_act"] += sum(max(s[1], 0.0) for s in segs)
        out["up_free"] += max(b.z - a.z, 0.0)

    for piece in rec["routes"]["F"]:
        load = sum(int(round(reqs[events[k][0]].q / STEP)) for k in piece if reqs[events[k][0]].kind == "D")
        pos = 0
        for k in piece:
            leg_add(pos, events[k][2], load * STEP)
            load += int(round(events[k][3] / STEP))
            pos = events[k][2]
        leg_add(pos, 0, load * STEP)
    e_act, e_free = out["act"].sum(), out["free"].sum()
    assert abs(e_act - rec["own"]["F"]) < 1e-6 * e_act, (path, e_act, rec["own"]["F"])
    cf = rec["cross"].get("F>FREE")
    if cf:
        assert abs(e_free - cf["energy"]) < 1e-6 * e_free, (path, e_free, cf["energy"])
    out["lam"], out["fam"], out["inst"] = lam, fam, rec["instance"]
    return out


def main():
    from joblib import Parallel, delayed
    files = sorted(glob.glob(os.path.join(CACHE, "study", "*.json")))
    res = Parallel(n_jobs=int(os.environ.get("STUDY_JOBS", "3")))(delayed(decompose_cell)(f) for f in files)
    L = ["# Engel etkisinin yatay/dikey ayrıştırması (KEŞİFSEL; decompose.py)", "",
         f"Her hücrenin saklı en iyi F rotası; her bacak gerçek yükte F'de en iyi seçeneğiyle ve aynı düğüm irtifalarında engelsiz düz uçuşla (FREE) kıyaslanır. {len(res)} hücre. "
         "H = yatay uçuş (P_h·ℓ/v), T = tırmanış/alçalış süresi için ek havada-kalma enerjisi, C = tırmanış (potansiyel) terimi m·g·Δz⁺/η_h. "
         "Fazlalık = gerçek − düz; bacak başına toplamsal. Kontroller: gerçek toplam = saklı F enerjisi; düz toplam = saklı F>FREE enerjisi (her hücrede sağlandı).", "",
         "| λ | F/FREE enerji artışı (medyan %) | ΔH payı | ΔT payı | ΔC payı | yatay yol uzaması (medyan %) | Σ tırmanış gerçek/düz (medyan) | bacak: doğrudan / üstten / çevreden (%) |",
         "|---|---|---|---|---|---|---|---|"]
    for lam in (1.0, 0.75, 0.5, 0.0):
        rs = [r for r in res if r["lam"] == lam]
        inc = [100 * (r["act"].sum() / r["free"].sum() - 1) for r in rs]
        ex = np.array([r["act"] - r["free"] for r in rs])
        tot = ex.sum(axis=1)
        share = ex / np.where(tot[:, None] == 0, np.nan, tot[:, None])
        det = [100 * (r["path_act"] / r["path_free"] - 1) for r in rs]
        upr = [r["up_act"] / r["up_free"] if r["up_free"] > 0 else np.nan for r in rs]
        kinds = {k: sum(r["kinds"][k] for r in rs) for k in ("direct", "over", "around", "free", "other")}
        nk = max(1, sum(kinds.values()))
        L.append(f"| {lam} | {np.median(inc):.1f} | {np.nanmedian(share[:, 0]):.2f} | {np.nanmedian(share[:, 1]):.2f} | {np.nanmedian(share[:, 2]):.2f} | "
                 f"{np.median(det):.1f} | {np.nanmedian(upr):.2f} | {100*kinds['direct']/nk:.0f} / {100*kinds['over']/nk:.0f} / {100*kinds['around']/nk:.0f} |")
    L += ["", "Not: paylar örnek başına hesaplanıp medyanlanır (toplamları 1 olmayabilir). ΔH<0 mümkündür (üstten geçiş yatay yolu kısaltabilir)."]
    open(os.path.join(RES, "DECOMP_auto.md"), "w", encoding="utf8").write(chr(10).join(L) + chr(10))
    print(chr(10).join(L))


if __name__ == "__main__":
    main()
