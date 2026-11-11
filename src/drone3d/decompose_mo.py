"""Horizontal / vertical decomposition of the obstacle effect under the Mo model (Round 10, finding E2; EXPLORATORY).

For the stored best Mo route of every cell (cache/study_mo/*.json) each leg (actual load) is evaluated with its best
geometric option under Mo, and compared with two obstacle-free counterfactuals at the SAME node altitudes:
  FREE  : one straight segment, climb and level flight simultaneous (t = max(...)), the model's own FREEMo leg;
  FREEq : the same straight leg flown as sequential V-H-V ((L,0) then (0,dz)), the SAME time structure as the 'around' option.
Excess(FREE)  = H + T + C terms difference; the FREEq comparison removes the V-H-V structural artefact.
Terms (Mo): H = P_h * ell / v ; T = P_h * (t - ell / v) ; C = max(p_c - P_h, 0) * up / v_up   (kappa_c = 1).
Also records how often the chosen 'around' option cruises at max(z_a, z_b) (street-level detour, f = 0).
Variant 'vdown5' (v_down = 5 m/s, UNSOURCED sensitivity, manual value is 3 m/s) re-chooses options for the same route.
Usage: python -I decompose_mo.py  -> 04_Results/DECOMP_Mo_auto.md
"""
import glob
import json
import math
import os
import sys
from dataclasses import replace

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from campaign import CACHE, RES, instance_specs  # noqa: E402
from energy import G, RHO, Hardware, segment_energy_wh  # noqa: E402
from exact import trip_events  # noqa: E402
from fast import STEP  # noqa: E402

V = 10.0


def terms(ell, dz, load_kg, hw):
    up, down = max(dz, 0.0), max(-dz, 0.0)
    m = hw.m0 + load_kg
    ph = hw.hover_power(load_kg)
    t = max(ell / V, up / hw.v_up, down / hw.v_down)
    thrust = m * G
    v_h = math.sqrt(thrust / (2.0 * RHO * hw.disk_area))
    p_c = thrust * (hw.v_up / 2.0 + math.sqrt((hw.v_up / 2.0) ** 2 + v_h ** 2)) / hw.eta_h
    return np.array([ph * (ell / V) / 3600.0, ph * (t - ell / V) / 3600.0, max(p_c - ph, 0.0) * up / hw.v_up / 3600.0])


def seg_total(segs, load_kg, hw):
    return sum(terms(e, d, load_kg, hw) for e, d in segs)


def decompose_cell(path, variant):
    from instances import DEPOT_A, make_instance
    from legs import LegCost
    from study import _depot_b, _env_region
    rec = json.load(open(path))
    fam, n, seed, lam = rec["instance"].split("_lam")[0].split("_")[0], rec["n"], rec["seed"], rec["lam"]
    region = "B" if fam == "RB" else "A"
    mix = (0.4, 0.3, 0.3)
    if fam == "M":
        mix = [s for s in instance_specs(["M"]) if s["n"] == 15 and s["seed"] == seed][0]["mix"]
    city, pluto, ras = _env_region(region)
    hw = Hardware() if variant == "default" else replace(Hardware(), v_down=5.0)
    inst = make_instance(city, pluto, DEPOT_A if region == "A" else _depot_b(city), n, seed, mix=mix, lam=lam)
    leg = LegCost(city, inst.nodes, model="Mo", raster=ras)
    reqs = inst.requests
    events = trip_events(reqs, tuple(range(len(reqs))))
    out = dict(act=np.zeros(3), free=np.zeros(3), freeq=np.zeros(3), kinds={}, around_f0=0, around_n=0,
               path_act=0.0, path_free=0.0, lam=lam, fam=fam)

    def leg_add(i, j, load_kg):
        if i == j:
            return
        best = None
        for kind, segs in leg.options(i, j):
            e = float(seg_total(segs, load_kg, hw).sum())
            if best is None or e < best[0]:
                best = (e, kind, segs)
        e, kind, segs = best
        a, b = inst.nodes[i], inst.nodes[j]
        L = math.hypot(b.x - a.x, b.y - a.y)
        out["act"] += seg_total(segs, load_kg, hw)
        out["free"] += terms(L, b.z - a.z, load_kg, hw)
        out["freeq"] += seg_total([(L, 0.0), (0.0, b.z - a.z)], load_kg, hw)
        out["kinds"][kind] = out["kinds"].get(kind, 0) + 1
        out["path_act"] += sum(s[0] for s in segs)
        out["path_free"] += L
        if kind == "around":
            out["around_n"] += 1
            h = a.z + segs[0][1]
            if abs(h - max(a.z, b.z)) < 1e-6:
                out["around_f0"] += 1

    for piece in rec["routes"]["Mo"]:
        load = sum(int(round(reqs[events[k][0]].q / STEP)) for k in piece if reqs[events[k][0]].kind == "D")
        pos = 0
        for k in piece:
            leg_add(pos, events[k][2], load * STEP)
            load += int(round(events[k][3] / STEP))
            pos = events[k][2]
        leg_add(pos, 0, load * STEP)
    if variant == "default":
        e_act = out["act"].sum()
        assert abs(e_act - rec["own"]["Mo"]) < 1e-6 * e_act, (path, e_act, rec["own"]["Mo"])
        cf = rec["cross"].get("Mo>FREEMo")
        if cf:
            assert abs(out["free"].sum() - cf["energy"]) < 1e-6 * out["free"].sum(), (path, out["free"].sum(), cf["energy"])
    return out


def summarize(res, label):
    L = [f"### Varyant: {label} ({len(res)} hücre)", "",
         "| λ | Mo/FREE artış (medyan %) | ΔH | ΔT | ΔC (FREE'ye göre pay, medyan) | Mo/FREEq artış (medyan %) | ΔH | ΔT | ΔC (FREEq'ya göre) | around bacaklarında f=0 payı | tür: doğrudan/üstten/çevreden (%) |",
         "|---|---|---|---|---|---|---|---|---|---|---|"]
    for lam in (1.0, 0.75, 0.5, 0.0):
        rs = [r for r in res if r["lam"] == lam]
        cells = []
        for key in ("free", "freeq"):
            inc = [100 * (r["act"].sum() / r[key].sum() - 1) for r in rs]
            ex = np.array([r["act"] - r[key] for r in rs])
            tot = ex.sum(axis=1)
            sh = ex / np.where(tot[:, None] == 0, np.nan, tot[:, None])
            cells += [f"{np.median(inc):.1f}"] + [f"{np.nanmedian(sh[:, k]):.2f}" for k in range(3)]
        f0 = sum(r["around_f0"] for r in rs) / max(1, sum(r["around_n"] for r in rs))
        kinds = {}
        for r in rs:
            for k, v in r["kinds"].items():
                kinds[k] = kinds.get(k, 0) + v
        nk = max(1, sum(kinds.values()))
        L.append(f"| {lam} | " + " | ".join(cells) + f" | {f0:.2f} | {100*kinds.get('direct', 0)/nk:.0f} / {100*kinds.get('over', 0)/nk:.0f} / {100*kinds.get('around', 0)/nk:.0f} |")
    return L


def main():
    from joblib import Parallel, delayed
    files = sorted(glob.glob(os.path.join(CACHE, os.environ.get("STUDY_IN","study_mo"), "*.json")))
    L = ["# Engel etkisinin yatay/dikey ayrıştırması, Mo modeli (KEŞİFSEL; decompose_mo.py; Round 10 E2)", "",
         "Her hücrenin saklı en iyi MO rotası; bacak başına en iyi seçenek ve aynı düğüm irtifalarında iki engelsiz karşıolgu: FREE (tek segment, eş zamanlı tırmanış) ve FREEq (ardışık V–H–V, 'around' ile aynı zaman yapısı). "
         "ΔH/ΔT/ΔC payları hücre başına hesaplanıp medyanlanır (toplamları 1 olmayabilir). v_down=5 yalnız KAYNAKSIZ duyarlılıktır (kılavuz değeri 3 m/s); κ_c=1 (seyir gücü=hover gücü) kaynaksız yapısal varsayım olarak sabit, değiştirilmedi.", ""]
    for variant, label in (("default", "ana (v_down=3)"), ("vdown5", "v_down=5 (kaynaksız duyarlılık)")):
        res = Parallel(n_jobs=int(os.environ.get("STUDY_JOBS", "8")))(delayed(decompose_cell)(f, variant) for f in files)
        L += summarize(res, label) + [""]
    open(os.path.join(RES, os.environ.get("DECOMP_OUT","DECOMP_Mo_auto.md")), "w", encoding="utf8").write(chr(10).join(L) + chr(10))


if __name__ == "__main__":
    main()


