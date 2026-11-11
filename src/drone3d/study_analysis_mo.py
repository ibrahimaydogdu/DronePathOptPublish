"""RQ2/RQ3 with the Mo model as MAIN model (Round 10 decision; EXPLORATORY labels as in STUDY_RQ123_auto).
Input: cache/study_mo/*.json (conditions Mo, 2D, NLMo, FREEMo).  Usage: python -I study_analysis_mo.py <dir> <out.md>
Regret = E_B(route_A) / E_B(reference) - 1; infeasible cross evaluations are excluded (counted).  Reference variants:
own = best route found under B itself; pooled = min over all routes evaluated in B's table.
"""
import glob
import json
import sys

import numpy as np
from scipy.stats import binomtest

LAMS = (1.0, 0.75, 0.5, 0.0)


def region(rec):
    f = rec["instance"].split("_")[0]
    return {"M": "A1", "RA": "A2", "RB": "B"}[f]


def regrets(recs, a, b, pooled):
    out = []
    for r in recs:
        c = r["cross"].get(f"{a}>{b}")
        if c is None:
            out.append(None)
            continue
        if not c["feasible"]:
            out.append(None)
            continue
        if pooled:
            cands = [r["cross"][f"{x}>{b}"]["energy"] for x in ("Mo", "2D", "NLMo", "FREEMo")
                     if f"{x}>{b}" in r["cross"] and r["cross"][f"{x}>{b}"]["feasible"]]
            ref = min(cands + [r["own"][b]])
        else:
            ref = r["own"][b]
        out.append(c["energy"] / ref - 1.0)
    return out


def main(d, out_md):
    recs = [json.load(open(f)) for f in sorted(glob.glob(d + "/*.json"))]
    noise = [v["rel_energy"] for r in recs for v in r["noise"].values()]
    tau = float(np.percentile(noise, 95))
    tau_mo = float(np.percentile([r["noise"]["Mo"]["rel_energy"] for r in recs if "Mo" in r["noise"]], 95))
    rng = np.random.default_rng(0)
    boot = [float(np.percentile(rng.choice(noise, len(noise)), 95)) for _ in range(2000)]
    tau_ci = (float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5)))
    L = ["# RQ2/RQ3, ANA MODEL = Mo (otomatik; study_analysis_mo.py; KEŞİFSEL)", "",
         f"Hücre sayısı: {len(recs)} (A1 {sum(region(r)=='A1' for r in recs)}, A2 {sum(region(r)=='A2' for r in recs)}, B {sum(region(r)=='B' for r in recs)}); K=8 koşu/koşul; "
         f"τ_E (gürültü tabanı %95; havuz = 4 koşul, bu veriden yeniden tahmin; ön-kayıtlı DEĞİL) = %{100*tau:.2f} (bootstrap %95 GA %{100*tau_ci[0]:.2f}–%{100*tau_ci[1]:.2f}); yalnız Mo koşulu: %{100*tau_mo:.2f}. 4/4 koşu bölmesi best-of-4 gürültüsüdür, pişmanlıkta kullanılan best-of-8 daha az gürültülü, yani τ_E muhafazakârdır. NLMo = Mo'nun yük-0 diliminin tüm yüklere kopyası; FREEMo = engelsiz düz uçuş (Mo).", "",
         "## Enerji düzeyi (kendi-en-iyi, medyan Wh) ve oranlar", "| λ | Mo | 2D | NLMo | FREEMo | Mo/2D | Mo/NLMo | Mo/FREEMo (%artış) |", "|---|---|---|---|---|---|---|---|"]
    for lam in LAMS:
        rs = [r for r in recs if r["lam"] == lam]
        m = lambda c: np.median([r["own"][c] for r in rs])
        L.append(f"| {lam} | {m('Mo'):.0f} | {m('2D'):.0f} | {m('NLMo'):.0f} | {m('FREEMo'):.0f} | "
                 f"{np.median([r['own']['Mo']/r['own']['2D'] for r in rs]):.2f} | {np.median([r['own']['Mo']/r['own']['NLMo'] for r in rs]):.3f} | "
                 f"{100*np.median([r['own']['Mo']/r['own']['FREEMo']-1 for r in rs]):.1f} |")
    for title, a in (("RQ2 — engel: FREEMo rotasının Mo'da pişmanlığı", "FREEMo"), ("RQ3 — yük: NLMo rotasının Mo'da pişmanlığı", "NLMo"),
                     ("2D rotasının Mo'da pişmanlığı", "2D")):
        L += ["", f"## {title}", "| λ | referans | n (uygun/toplam) | medyan % | pişmanlık>τ_E payı | işaret testi (pozitif/n, p; yalnız kendi en iyi referans) | uygunsuz dahil üst sınır (pay) | sefer sayısı Mo / diğer |", "|---|---|---|---|---|---|---|---|"]
        for lam in LAMS:
            rs = [r for r in recs if r["lam"] == lam]
            for pooled in (False, True):
                g = regrets(rs, a, "Mo", pooled)
                v = np.array([x for x in g if x is not None])
                pos = int((v > 0).sum())
                p = binomtest(pos, len(v), 0.5).pvalue if len(v) else float("nan")
                ptxt = "—" if pooled else ("p<1e-4" if p < 1e-4 else f"p={p:.4f}")
                ub = (np.sum(v > tau) + sum(x is None for x in g)) / len(g)
                ta = np.mean([r["cross"][f"{a}>Mo"]["trips_b"] for r in rs if f"{a}>Mo" in r["cross"]])
                tb = np.mean([r["cross"][f"{a}>Mo"]["trips_a"] for r in rs if f"{a}>Mo" in r["cross"]])
                L.append(f"| {lam} | {'toplanmış' if pooled else 'kendi en iyi'} | {len(v)}/{len(g)} | {100*np.median(v):.1f} | {np.mean(v > tau):.2f} | {'—' if pooled else str(pos)+'/'+str(len(v))+', '+ptxt} | {ub:.2f} | {ta:.2f} / {tb:.2f} |")
    L += ["", "## Bölgeye göre (toplanmış referans, NLMo→Mo ve FREEMo→Mo medyan % / τ_E aşan pay)", "| bölge | λ | NLMo→Mo | FREEMo→Mo | 2D→Mo |", "|---|---|---|---|---|"]
    for reg in ("A1", "A2", "B"):
        for lam in LAMS:
            rs = [r for r in recs if r["lam"] == lam and region(r) == reg]
            cells = []
            for a in ("NLMo", "FREEMo", "2D"):
                v = np.array([x for x in regrets(rs, a, "Mo", True) if x is not None])
                cells.append(f"{100*np.median(v):.1f} / {np.mean(v > tau):.2f} (n={len(v)})")
            L.append(f"| {reg} | {lam} | " + " | ".join(cells) + " |")

    noiseY = [v["d_struct"] for r in recs for v in r["noise"].values()]
    tauY = float(np.percentile(noiseY, 95))
    L += ["", f"## Ön-kayıtlı kural (Mo): pişmanlık > τ_E VE d_struct > τ_Y; τ_E=%{100*tau:.2f}, τ_Y={tauY:.3f}", "| λ | çift | n | pişmanlık>τ_E | d_struct>τ_Y | ikisi birden |", "|---|---|---|---|---|---|"]
    for lam in LAMS:
        rs = [r for r in recs if r["lam"] == lam]
        for a in ("FREEMo", "NLMo", "2D"):
            rg, ds = [], []
            for r, g in zip(rs, regrets(rs, a, "Mo", True)):
                if g is None:
                    continue
                rg.append(g); ds.append(r["cross"][f"{a}>Mo"]["d_struct"])
            rg, ds = np.array(rg), np.array(ds)
            L.append(f"| {lam} | {a}→Mo | {len(rg)} | {np.mean(rg>tau):.2f} | {np.mean(ds>tauY):.2f} | {np.mean((rg>tau)&(ds>tauY)):.2f} |")

    L += ["", "## Birincil işaret testi: örnek başına (28 bağımsız birim), λ-ortalamalı, kendi en iyi referans, uygun λ'lar", "| çift | n örnek | pozitif | medyan % | en küçük % | p (iki taraflı binom) |", "|---|---|---|---|---|---|"]
    for a in ("FREEMo", "NLMo", "2D"):
        per = {}
        for r in recs:
            c = r["cross"].get(f"{a}>Mo")
            if c is None or not c["feasible"]:
                continue
            per.setdefault(r["instance"].split("_lam")[0], []).append(c["energy"] / r["own"]["Mo"] - 1.0)
        m = np.array([np.mean(v) for v in per.values()])
        pos = int((m > 0).sum()); p = binomtest(pos, len(m), 0.5).pvalue
        L.append(f"| {a}→Mo | {len(m)} | {pos} | {100*np.median(m):.1f} | {100*m.min():.1f} | {'<1e-4' if p < 1e-4 else f'{p:.4f}'} |")
    L += ["", "Not (E1, Round 10): d_struct eşiği τ_Y rastgele 3-sefer bölmesinin beklenen d_struct'ından (≈0,44) büyük; yapı kolu bu çalışmada SINAMA GÜCÜNE SAHİP DEĞİL, 'ikisi birden' payı fiilen d_struct payına eşittir. Sonuç 'ayırt edilemez' değil 'sınanamadı'. Mo ve 2D sütunları eski F koşusuyla aynı rotalardır (yeni veri değil); yalnız NLMo/FREEMo yenidir. FREEMo/NLMo kendi-en-iyi enerjileri Mo'nun alt sınırıdır (serbest dikey/yüksüz varsayım)."]
    L += ["", "Not: negatif pişmanlık çözücü gürültüsüdür; uygunsuz çapraz değerlendirmeler dışlanır (AŞAĞI yanlı). Bölgelerde ve λ'larda örnekler bağımsız değildir; işaret testi düzeltmesizdir."]
    open(out_md, "w", encoding="utf8").write(chr(10).join(L) + chr(10))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
