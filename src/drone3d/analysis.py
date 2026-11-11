"""Campaign analysis (ANALYSIS_PLAN.md v2). Pure functions + a report writer.

Usage: python -I analysis.py <runs.csv> <instances.csv> <out.md>
"""
from __future__ import annotations
import itertools
import sys

import numpy as np

MAIN = ("TLBO", "GWO", "ABC")
BASE = ("RS", "HC", "NN2OPT")
DELTA_EQ = 0.02


def a12(x, y) -> float:
    """P(x < y) + 0.5 P(x == y): probability that x (lower energy) beats y."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    gt = (x[:, None] < y[None, :]).sum()
    eq = (x[:, None] == y[None, :]).sum()
    return float((gt + 0.5 * eq) / (len(x) * len(y)))


def holm(pvals: list[float]) -> list[float]:
    m = len(pvals)
    order = np.argsort(pvals)
    adj = np.empty(m)
    run = 0.0
    for rank, idx in enumerate(order):
        run = max(run, (m - rank) * pvals[idx])
        adj[idx] = min(1.0, run)
    return adj.tolist()


def rel_iqr(v) -> float:
    v = np.asarray(v, float)
    med = np.median(v)
    return float((np.percentile(v, 75) - np.percentile(v, 25)) / med) if med > 0 else 0.0


def struggle(ml: list[dict], s: list[dict],
             delta=0.02, iqr_mult=2.0, a12_min=0.64, frac_unreliable=0.5, frac_value=0.5,
             s_hit=0.8, s_frac=0.25) -> dict:
    """ANALYSIS_PLAN section 6 as ONE Boolean expression.

    ml: per M/L instance dict(medians={alg: median energy}, q=rel IQR of best algorithm, b=A12 of best alg vs best baseline)
    s : per S instance dict(hit=max over algorithms of optimum-hit rate)
    """
    unreliable, adds = [], []
    for inst in ml:
        med = np.array(list(inst["medians"].values()), float)
        r = (med.max() - med.min()) / med.min()
        unreliable.append(r > max(delta, iqr_mult * inst["q"]))
        adds.append(inst["b"] >= a12_min)
    struggle_ml = bool(ml) and (np.mean(unreliable) > frac_unreliable) and (np.mean(adds) >= frac_value)
    struggle_s = bool(s) and (np.mean([inst["hit"] < s_hit for inst in s]) > s_frac)
    return {"STRUGGLE_ML": bool(struggle_ml), "STRUGGLE_S": bool(struggle_s),
            "STRUGGLE": bool(struggle_ml or struggle_s),
            "mean_unreliable": float(np.mean(unreliable)) if ml else None,
            "mean_adds_value": float(np.mean(adds)) if ml else None}


def struggle_inputs(runs):
    """Per-instance inputs for struggle() computed from the raw run table (ANALYSIS_PLAN section 6).

    runs: DataFrame with columns instance, family, algorithm, feasible, energy, exact (S only).
    ML families M and L; S families S1 and S2. NN2OPT is a single deterministic value: its A12 is the share of
    a*'s runs below that value (ties counted half); no hypothesis test is run for it.
    """
    ml, s = [], []
    for inst, g in runs.groupby("instance"):
        fam = g["family"].iloc[0]
        if fam in ("M", "L"):
            ok = g[g["feasible"] == True]  # noqa: E712
            meds = {a: float(np.median(ok[ok["algorithm"] == a]["energy"])) for a in MAIN
                    if (ok["algorithm"] == a).any()}
            if len(meds) < len(MAIN):
                continue
            best_alg = min(meds, key=meds.get)
            x = ok[ok["algorithm"] == best_alg]["energy"].values
            cands = {}
            for b in BASE:
                v = ok[ok["algorithm"] == b]["energy"].values
                if len(v):
                    cands[b] = v
            if not cands:
                continue
            best_base = min(cands, key=lambda b: float(np.median(cands[b])))
            ml.append({"instance": inst, "medians": meds, "q": rel_iqr(x),
                       "b": a12(x, cands[best_base]), "best_alg": best_alg, "best_base": best_base})
        elif fam in ("S1", "S2"):
            ex = g["exact"].iloc[0]
            if ex == "" or pd_isna(ex):
                continue
            ex = float(ex)
            hits = []
            for a in MAIN:
                ga = g[g["algorithm"] == a]
                hits.append(float(np.mean([(r.feasible and r.energy <= ex * (1 + 1e-6)) for r in ga.itertuples()])) if len(ga) else 0.0)
            s.append({"instance": inst, "hit": max(hits)})
    return ml, s


def pd_isna(v):
    try:
        return v != v
    except Exception:
        return False


def mann_whitney(x, y):
    from scipy.stats import mannwhitneyu
    if np.all(np.asarray(x) == np.asarray(x)[0]) and np.all(np.asarray(y) == np.asarray(x)[0]):
        return 1.0
    return float(mannwhitneyu(x, y, alternative="two-sided").pvalue)


def family_tests(med: dict[str, np.ndarray]):
    """med[alg] = per-instance medians (aligned). Friedman then Holm-Wilcoxon on the 3 pairs."""
    from scipy.stats import friedmanchisquare, wilcoxon
    algs = list(med)
    out = {"friedman_p": float(friedmanchisquare(*[med[a] for a in algs]).pvalue), "pairs": []}
    raw = []
    for a, b in itertools.combinations(algs, 2):
        d = np.asarray(med[a]) - np.asarray(med[b])
        d = d[d != 0]
        raw.append(1.0 if len(d) == 0 else float(wilcoxon(d).pvalue))
    adj = holm(raw)
    for (a, b), p, q in zip(itertools.combinations(algs, 2), raw, adj):
        out["pairs"].append((a, b, p, q))
    return out


def main(runs_csv, inst_csv, out_md, long_csv=None):
    """Round 06 fixes: medians only over the instance set common to the compared algorithms (NN2OPT reported by
    feasibility count, never by a median over its feasible subset); best-known includes the long runs; S hit rate by D;
    STRUGGLE wired in."""
    import os
    import pandas as pd
    runs = pd.read_csv(runs_csv)
    insts = pd.read_csv(inst_csv)
    lines = ["# Kampanya analizi (otomatik üretildi; ANALYSIS_PLAN.md v2 + Round 06 düzeltmeleri)", ""]
    long = pd.read_csv(long_csv) if long_csv and os.path.exists(long_csv) else None
    for fam, g in runs.groupby("family"):
        lines.append(f"## Aile {fam}")
        feas = g.groupby("algorithm")["feasible"].mean()
        lines.append("Uygun koşu oranı: " + ", ".join(f"{a} {v:.3f}" for a, v in feas.items()))
        nn = g[g["algorithm"] == "NN2OPT"]
        if len(nn):
            lines.append(f"NN+2-opt uygun örnek: {int(nn['feasible'].sum())}/{len(nn)} (medyan özetine GİRMEZ; seçilim yanlılığı)")
        gg = g[(g["feasible"] == True) & (g["algorithm"] != "NN2OPT")]  # noqa: E712
        med = gg.groupby(["instance", "algorithm"])["energy"].median().unstack()
        med = med.dropna(subset=[a for a in MAIN if a in med.columns])
        lines.append(f"Ortak örnek sayısı: {len(med)}. Medyan enerji özeti (Wh): " +
                     "; ".join(f"{a}: {med[a].median():.1f}" for a in med.columns))
        if all(a in med.columns for a in MAIN) and len(med) >= 3:
            ft = family_tests({a: med[a].values for a in MAIN})
            lines.append(f"Friedman p = {ft['friedman_p']:.4g}; " + "; ".join(
                f"{a}-{b}: p={p:.4g}, Holm={q:.4g}" for a, b, p, q in ft["pairs"]))
            lines.append("Not: Wilcoxon p değeri örnek sayısının verebileceği en küçük değerdir; kanıt 'her örnekte aynı yön'dür.")
        if fam in ("M", "L") and long is not None:
            lg = long[(long["family"] == fam) & (long["feasible"] == True)]  # noqa: E712
            ratios = []
            for inst, gi in g[g["feasible"] == True].groupby("instance"):  # noqa: E712
                li = lg[lg["instance"] == inst]
                if len(li):
                    ratios.append(1 - li["energy"].min() / gi["energy"].min())
            if ratios:
                lines.append(f"En iyi bilinen (uzun koşu dahil) 1× en iyiden medyan %{100*np.median(ratios):.1f} düşük "
                             f"({len(ratios)} örnek): M/L'de 'optimuma boşluk' iddia edilmez.")
        lines.append("")
    s_runs = runs[runs["family"].isin(("S1", "S2"))].merge(insts[["instance", "events"]], on="instance", how="left")
    if len(s_runs):
        s_runs["hit"] = s_runs["feasible"] & (s_runs["energy"] <= s_runs["exact"] * (1 + 1e-6))
        t = s_runs[s_runs["algorithm"].isin(MAIN)].groupby(["events", "algorithm"])["hit"].mean().unstack().round(2)
        lines += ["## S ailesi: optimuma ulaşma oranı, D (olay sayısı) bazında", t.to_string(), ""]
    ml, s = struggle_inputs(runs)
    st = struggle(ml, s)
    lines += ["## Zorlanma göstergesi (ANALYSIS_PLAN §6)", str(st),
              "Yorum: tetikleyici araç uyuşmazlığını ölçer; 'problem zor' kanıtı değildir (Round 06).", ""]
    open(out_md, "w", encoding="utf8").write(chr(10).join(lines))


if __name__ == "__main__":
    main(*sys.argv[1:5])
