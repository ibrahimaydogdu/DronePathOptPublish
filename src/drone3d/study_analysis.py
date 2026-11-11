"""Analysis of study.py output per PARAMETRIC_PLAN section 6 (rules fixed before data).
Usage: python -I study_analysis.py <cache/study dir> <out.md>"""
import glob
import os
sys_path_dir = os.path.dirname(os.path.abspath(__file__))
import json
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

MODELS = ("F", "D", "Mo", "L", "S0", "2D")


def main(d, out):
    recs = [json.load(open(f)) for f in sorted(glob.glob(d + "/*.json"))]
    noiseE = [v["rel_energy"] for r in recs for v in r["noise"].values()]
    noiseY = [v["d_struct"] for r in recs for v in r["noise"].values()]
    tauE, tauY = float(np.percentile(noiseE, 95)), float(np.percentile(noiseY, 95))
    L = ["# Model-biçimi çalışması: çıktı (otomatik; PARAMETRIC_PLAN §6)", "",
         f"Hücre sayısı (örnek×λ): {len(recs)}; koşu sayısı/koşul: K=8. Gürültü tabanı (%95): τ_E = {100*tauE:.2f} %, τ_Y = {tauY:.3f} (d_struct). "
         f"Aynı-model yarı-bölme medyanları: enerji farkı {100*np.median(noiseE):.2f} %, d_struct {np.median(noiseY):.3f}.", ""]
    lams = sorted({r["lam"] for r in recs}, reverse=True)

    def pair_stats(a, b, lam=None):
        rs = [r for r in recs if (lam is None or r["lam"] == lam) and f"{a}>{b}" in r["cross"]
              and r["cross"][f"{a}>{b}"]["feasible"]]      # Round 07: infeasible cross-evaluations are excluded, counted separately
        reg = np.array([r["cross"][f"{a}>{b}"]["regret"] for r in rs])
        ds = np.array([r["cross"][f"{a}>{b}"]["d_struct"] for r in rs])
        n_all = sum(1 for r in recs if (lam is None or r["lam"] == lam) and f"{a}>{b}" in r["cross"])
        inf = np.array([n_all - len(rs)] * 1) / max(n_all, 1)
        sens = (reg > tauE) & (ds > tauY)
        return dict(n=len(rs), reg_med=float(np.median(reg)) if len(rs) else float("nan"),
                    reg_max=float(reg.max()) if len(rs) else float("nan"), d_med=float(np.median(ds)) if len(rs) else float("nan"),
                    sens=float(sens.mean()) if len(rs) else float("nan"), infeas=float(inf.mean()))

    L.append("## RQ1 — çapraz pişmanlık (A'nın rotası B'de): λ başına, satır A → sütun B, hücre = medyan pişmanlık % / 'duyarlı' örnek payı")
    for lam in lams:
        L += ["", f"### λ = {lam}", "| A \\ B | " + " | ".join(MODELS) + " |", "|---|" + "---|" * len(MODELS)]
        for a in MODELS:
            row = []
            for b in MODELS:
                if a == b:
                    row.append("–")
                else:
                    s = pair_stats(a, b, lam)
                    row.append(f"{100*s['reg_med']:.1f} / {s['sens']:.2f}")
            L.append(f"| {a} | " + " | ".join(row) + " |")
    L += ["", "## Model düzeyinde özet (tüm λ, A≠B): duyarlı çift payı (pişmanlık>τ_E VE d_struct>τ_Y)"]
    for lam in lams:
        vals = [pair_stats(a, b, lam)["sens"] for a in MODELS for b in MODELS if a != b]
        L.append(f"- λ={lam}: ortalama duyarlı pay {np.mean(vals):.2f}; en büyük medyan pişmanlık "
                 f"{100*max(pair_stats(a, b, lam)['reg_med'] for a in MODELS for b in MODELS if a != b):.1f} %")
    L += ["", "## Enerji düzeyi etkisi (rota aynı olsa da enerji farkı): λ başına, kendi-en-iyi enerjilerin medyanı (Wh)",
          "| λ | " + " | ".join(MODELS + ("NL", "FREE")) + " |", "|---|" + "---|" * (len(MODELS) + 2)]
    for lam in lams:
        rs = [r for r in recs if r["lam"] == lam]
        L.append(f"| {lam} | " + " | ".join(f"{np.median([r['own'][c] for r in rs]):.0f}" for c in MODELS + ("NL", "FREE")) + " |")

    L += ["", "## RQ2 — bina engeli (F, engelli) ve engelsiz düz uçuş (FREE)", "| λ | enerji artışı engelden (medyan %) | FREE rotası F'de pişmanlık (medyan %) | F rotası FREE'de pişmanlık | duyarlı pay (FREE→F) | F rotasının sefer sayısı / FREE |", "|---|---|---|---|---|---|"]
    for lam in lams:
        rs = [r for r in recs if r["lam"] == lam]
        up = np.median([100 * (r["own"]["F"] / r["own"]["FREE"] - 1) for r in rs])
        s1, s2 = pair_stats("FREE", "F", lam), pair_stats("F", "FREE", lam)
        ta = np.mean([r["cross"]["F>FREE"]["trips_a"] for r in rs]); tb = np.mean([r["cross"]["F>FREE"]["trips_b"] for r in rs])
        L.append(f"| {lam} | {up:.1f} | {100*s1['reg_med']:.1f} | {100*s2['reg_med']:.1f} | {s1['sens']:.2f} | {ta:.2f} / {tb:.2f} |")

    L += ["", "## RQ3 — yük etkisi (NL: yükten bağımsız enerji; F: yüke bağlı)", "| λ | F enerjisi / NL enerjisi (medyan oran) | NL rotası F'de pişmanlık (medyan %) | duyarlı pay (NL→F) | sefer sayısı F / NL |", "|---|---|---|---|---|"]
    for lam in lams:
        rs = [r for r in recs if r["lam"] == lam]
        ratio = np.median([r["own"]["F"] / r["own"]["NL"] for r in rs])
        s = pair_stats("NL", "F", lam)
        ta = np.mean([r["cross"]["F>NL"]["trips_a"] for r in rs]); tb = np.mean([r["cross"]["F>NL"]["trips_b"] for r in rs])
        L.append(f"| {lam} | {ratio:.3f} | {100*s['reg_med']:.1f} | {s['sens']:.2f} | {ta:.2f} / {tb:.2f} |")

    L += exploratory(recs, tauE)
    L += pooled(recs, tauE)
    L += ["", "## Çözücü güvenilirliği", f"Uygun koşu sayısı (hücre×koşul başına 8): min {min(min(r['feasible_runs'].values()) for r in recs)}.",
          "Batarya/kapasite ihlali (A rotası B'de uygun değil) payı, yalnız A≠B çiftleri: "
          f"{100*np.mean([not v['feasible'] for r in recs for k, v in r['cross'].items() if k.split('>')[0] != k.split('>')[1]]):.1f} % "
          "(pişmanlık istatistikleri uygunsuz çapraz değerlendirmeleri DIŞLAR)"]
    open(out, "w", encoding="utf8").write(chr(10).join(L) + chr(10))

def pooled(recs, tauE):
    """KESIFSEL (Round 07 next step): pooled reference. For condition B the best known energy is the minimum over all
    routes (own and every other condition's) evaluated in B's table, which removes the negative-regret artefact."""
    L = ["", "## KEŞİFSEL: toplanmış referansla pişmanlık (referans = B tablosunda tüm rotaların en iyisi)", "",
         "| bölge | λ | RQ1 (5 yük-bağımlı model arası, A≠B) medyan % / τ_E'yi aşan pay | NL→F medyan % / pay | 2B→F medyan % / pay | FREE→F medyan % / pay |", "|---|---|---|---|---|---|"]
    regions = {}
    for r in recs:
        reg = "A2" if r["instance"].startswith("RA") else ("B" if r["instance"].startswith("RB") else "A1")
        regions.setdefault(reg, []).append(r)
    regions["tümü"] = recs
    def pr(rs, a, b):
        out = []
        for r in rs:
            cands = [r["cross"][f"{x}>{b}"]["energy"] for x in ("F", "D", "Mo", "L", "S0", "2D", "NL", "FREE")
                     if f"{x}>{b}" in r["cross"] and r["cross"][f"{x}>{b}"]["feasible"]]
            c = r["cross"].get(f"{a}>{b}")
            if cands and c and c["feasible"]:
                out.append(c["energy"] / min(cands) - 1.0)
        return np.array(out)
    for reg, rs in regions.items():
        for lam in sorted({r["lam"] for r in rs}, reverse=True):
            q = [r for r in rs if r["lam"] == lam]
            core = np.concatenate([pr(q, a, b) for a in ("F", "D", "Mo", "L", "S0") for b in ("F", "D", "Mo", "L", "S0") if a != b])
            def cell(a, b):
                v = pr(q, a, b)
                return f"{100*np.median(v):.1f} / {np.mean(v > tauE):.2f} (n={len(v)})" if len(v) else "–"
            L.append(f"| {reg} | {lam} | {100*np.median(core):.1f} / {np.mean(core > tauE):.2f} | {cell('NL', 'F')} | {cell('2D', 'F')} | {cell('FREE', 'F')} |")
    return L


def _request_labels(rec, events):
    lab = {}
    return lab


def ari(la, lb):
    """Adjusted Rand index between two trip labelings (0 = chance, 1 = identical)."""
    from math import comb
    n = len(la)
    ca, cb, cab = {}, {}, {}
    for x, y in zip(la, lb):
        ca[x] = ca.get(x, 0) + 1; cb[y] = cb.get(y, 0) + 1; cab[(x, y)] = cab.get((x, y), 0) + 1
    sa = sum(comb(v, 2) for v in ca.values()); sb = sum(comb(v, 2) for v in cb.values())
    sab = sum(comb(v, 2) for v in cab.values())
    tot = comb(n, 2)
    exp = sa * sb / tot if tot else 0.0
    mx = 0.5 * (sa + sb)
    return 1.0 if mx == exp else (sab - exp) / (mx - exp)


def exploratory(recs, tauE):
    """KESIFSEL (Round 07; added after seeing results): components of the decision rule, chance-corrected structure agreement."""
    import os
    from campaign import _env
    from exact import trip_events
    from instances import DEPOT_A, make_instance
    city, pluto, _ = _env()
    reqs_by_seed = {}
    for r in recs:
        if r["seed"] not in reqs_by_seed:
            reqs_by_seed[r["seed"]] = make_instance(city, pluto, DEPOT_A, r["n"], r["seed"], mix=(0.4, 0.3, 0.3)).requests
    L = ["", "## KEŞİFSEL (Round 07 sonrası eklendi; kural bileşenleri ayrı ayrı)"]
    L += ["", "### RQ3 bileşenleri (NL rotası F tablosunda; yalnız uygun çapraz değerlendirmeler)",
          "| λ | hücre (uygun/8) | pişmanlık>τ_E payı | pişmanlık>2τ_E payı | işaret testi (pozitif/n) | F sefer sayısı − NL (ort.) |", "|---|---|---|---|---|---|"]
    from scipy.stats import binomtest
    for lam in sorted({r["lam"] for r in recs}, reverse=True):
        rs = [r for r in recs if r["lam"] == lam]
        ok = [r for r in rs if r["cross"]["NL>F"]["feasible"]]
        reg = np.array([r["cross"]["NL>F"]["regret"] for r in ok])
        pos = int((reg > 0).sum())
        p = binomtest(pos, len(reg), 0.5).pvalue if len(reg) else float("nan")
        dt = np.mean([r["cross"]["F>NL"]["trips_a"] - r["cross"]["F>NL"]["trips_b"] for r in rs])
        L.append(f"| {lam} | {len(ok)}/{len(rs)} | {np.mean(reg > tauE):.2f} | {np.mean(reg > 2*tauE):.2f} | {pos}/{len(reg)} (p={p:.3f}) | {dt:+.2f} |")
    L += ["", "### Şansa düzeltilmiş sefer-yapısı uyumu (ARI, 0 = rastgele, 1 = aynı; en iyi rotalar arası; medyan)",
          "| λ | F–D | F–Mo | F–L | F–S0 | F–2D | F–NL | F–FREE |", "|---|---|---|---|---|---|---|---|"]
    for lam in sorted({r["lam"] for r in recs}, reverse=True):
        row = []
        for b in ("D", "Mo", "L", "S0", "2D", "NL", "FREE"):
            v = []
            for r in [x for x in recs if x["lam"] == lam]:
                reqs = reqs_by_seed[r["seed"]]
                ev = trip_events(reqs, tuple(range(len(reqs))))
                def lab(pieces):
                    t = {}
                    for ti, pc in enumerate(pieces):
                        for k in pc:
                            t[ev[k][0]] = ti
                    return [t[i] for i in range(len(reqs))]
                if r["routes"]["F"] and r["routes"][b]:
                    v.append(ari(lab(r["routes"]["F"]), lab(r["routes"][b])))
            row.append(f"{np.median(v):.2f}")
        L.append(f"| {lam} | " + " | ".join(row) + " |")
    L += ["", "### Enerji düzeyi: F/2D oranı (örnek başına) aralığı, λ başına",
          "| λ | min | medyan | maks |", "|---|---|---|---|"]
    for lam in sorted({r["lam"] for r in recs}, reverse=True):
        q = [r["own"]["F"] / r["own"]["2D"] for r in recs if r["lam"] == lam]
        L.append(f"| {lam} | {min(q):.2f} | {np.median(q):.2f} | {max(q):.2f} |")
    L += ["", "Not: d_order raporlanmıyor (bilgi taşımıyor, Round 07). Negatif pişmanlık çözücü gürültüsüdür; bütün sonuçlar K=8 en iyi-koşu çözünürlüğüyle (≈%1–2) sınırlıdır."]
    return L


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
