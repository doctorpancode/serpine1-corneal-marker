# -*- coding: utf-8 -*-
"""Revise-support analyses for 草稿1.2 (biomarker reframing).
Produces verified numbers to replace the unverified ones in 1.1.
"""
import json, math, re
import numpy as np
import pandas as pd
from itertools import combinations

RNG = np.random.default_rng(20261111)
OUT = {}

# ---------------------------------------------------------------- 1. GSE241715
p = r"D1执行结果\第一天执行结果\data\processed\GSE241715_expression.tsv"
ex = pd.read_csv(p, sep="\t", index_col=0)
ss = pd.read_csv(r"D1执行结果\第一天执行结果\data\processed\GSE241715_sample_sheet.tsv", sep="\t")
samples = [c for c in ex.columns]
grp = dict(zip(ss["sample"], ss["group"]))
genes = ["SERPINE1", "TMSB4X", "P4HB", "GJA1", "S100A13", "SOD1", "PRDX2", "VDAC1", "CAPG", "XIAP"]
genes = [g for g in genes if g in ex.index]

def auc(x_pos, x_neg):
    """Mann-Whitney AUC."""
    x_pos = np.asarray(x_pos, float); x_neg = np.asarray(x_neg, float)
    gt = sum((a > b) for a in x_pos for b in x_neg)
    eq = sum((a == b) for a in x_pos for b in x_neg)
    return (gt + 0.5 * eq) / (len(x_pos) * len(x_neg))

def boot_auc(pos, neg, n=4000):
    pos = np.asarray(pos, float); neg = np.asarray(neg, float)
    vals = []
    for _ in range(n):
        vals.append(auc(RNG.choice(pos, len(pos), True), RNG.choice(neg, len(neg), True)))
    return np.percentile(vals, [2.5, 97.5])

stats241 = {}
for g in genes:
    v = ex.loc[g]
    d = {}
    for gname in ["Control", "Keratoconus", "Viral", "Bacterial_fungal"]:
        cols = [s for s in samples if grp.get(s) == gname]
        x = v[cols].astype(float).values
        d[gname] = dict(mean=float(np.mean(x)), sd=float(np.std(x, ddof=1)), n=len(x))
    ctrl = [s for s in samples if grp.get(s) == "Control"]
    les = [s for s in samples if grp.get(s) != "Control"]
    bf = [s for s in samples if grp.get(s) == "Bacterial_fungal"]
    d["auc_lesion"] = float(auc(v[les].values, v[ctrl].values))
    d["auc_lesion_ci"] = [float(x) for x in boot_auc(v[les].values, v[ctrl].values)]
    d["auc_bf"] = float(auc(v[bf].values, v[ctrl].values))
    d["auc_bf_ci"] = [float(x) for x in boot_auc(v[bf].values, v[ctrl].values)]
    d["fold_ctrl_to_bf"] = float(d["Bacterial_fungal"]["mean"] / d["Control"]["mean"]) if d["Control"]["mean"] else None
    stats241[g] = d
OUT["GSE241715"] = stats241

# LOOCV logistic (single predictor, unpenalised -> guarded)
def loocv_auc(x, y):
    x = np.asarray(x, float); y = np.asarray(y, int)
    pred = np.zeros(len(x))
    ok = True
    for i in range(len(x)):
        m = np.ones(len(x), bool); m[i] = False
        xt, yt = x[m], y[m]
        # simple logistic by IRLS
        b = np.zeros(2)
        X = np.column_stack([np.ones(len(xt)), xt])
        for _ in range(60):
            eta = X @ b
            pr = 1 / (1 + np.exp(-eta))
            W = pr * (1 - pr) + 1e-10
            try:
                b_new = b + np.linalg.solve(X.T @ (X * W[:, None]), X.T @ (yt - pr))
            except np.linalg.LinAlgError:
                ok = False; break
            if np.max(np.abs(b_new - b)) < 1e-8:
                b = b_new; break
            b = b_new
        pred[i] = 1 / (1 + np.exp(-(b[0] + b[1] * x[i])))
    if not ok:
        return None, False
    yv = y
    return float(auc(pred[yv == 1], pred[yv == 0])), True

ctrl_cols = [s for s in samples if grp.get(s) == "Control"]
les_cols = [s for s in samples if grp.get(s) != "Control"]
y = np.array([0] * len(ctrl_cols) + [1] * len(les_cols))
for g in ["SERPINE1", "TMSB4X"]:
    x = np.concatenate([ex.loc[g, ctrl_cols].astype(float).values, ex.loc[g, les_cols].astype(float).values])
    a, ok = loocv_auc(x, y)
    OUT.setdefault("LOOCV", {})[g] = dict(auc=a, ok=ok)

# ---------------------------------------------------------------- 2. GSE58291
raw = pd.read_csv(r"第9天执行结果\data\GSE58291_expr.tsv", sep="\t", low_memory=False)
avg_cols = [c for c in raw.columns if c.endswith(".AVG_Signal")]
smap = {}
for c in avg_cols:
    nm = re.sub(r"\s+", "", c.replace(".AVG_Signal", ""))  # "M 23" -> "M23"
    smap[c] = nm.strip()
raw["SYM"] = raw.iloc[:, 1].astype(str).str.strip().str.upper()
sub = raw[raw["SYM"].isin(["SERPINE1", "TMSB4X", "P4HB", "GJA1"])]
sub = sub.groupby("SYM")[avg_cols].max()
sgrp = pd.read_csv(r"第9天执行结果\results\D9_GSE58291_sample_group.tsv", sep="\t")
gmap = dict(zip(sgrp["sample"], sgrp["group"]))
res58291 = {}
for g in sub.index:
    vals = {smap[c]: float(sub.loc[g, c]) for c in avg_cols if smap[c] in gmap}
    ctr = [v for k, v in vals.items() if gmap[k] == "Control"]
    bk = [v for k, v in vals.items() if gmap[k] == "BK"]
    fk = [v for k, v in vals.items() if gmap[k] == "FK"]
    res58291[g] = dict(
        n=dict(Control=len(ctr), FK=len(fk), BK=len(bk)),
        mean=dict(Control=float(np.mean(ctr)), FK=float(np.mean(fk)), BK=float(np.mean(bk))),
        median=dict(Control=float(np.median(ctr)), FK=float(np.median(fk)), BK=float(np.median(bk))),
        auc_bk_vs_ctrl=float(auc(bk, ctr)),
        auc_bk_vs_ctrl_ci=[float(x) for x in boot_auc(bk, ctr)],
        auc_inf_vs_ctrl=float(auc(bk + fk, ctr)),
        auc_inf_vs_ctrl_ci=[float(x) for x in boot_auc(bk + fk, ctr)],
        delta_bk=float(np.mean(bk) - np.mean(ctr)),
        delta_fk=float(np.mean(fk) - np.mean(ctr)),
    )
OUT["GSE58291"] = res58291

# ---------------------------------------------------------------- 3. DEG counts
deg = pd.read_csv(r"第3天执行结果\results\D3_limma_deg.tsv", sep="\t")
deg.columns = [c.strip() for c in deg.columns]
deg = deg.dropna(subset=["adjP"])
OUT["DEG"] = dict(
    total=int(len(deg)),
    padj05=int((deg["adjP"] < 0.05).sum()),
    strict=int(((deg["adjP"] < 0.05) & (deg["logFC"].abs() > 0.5)).sum()),
    up_strict=int(((deg["adjP"] < 0.05) & (deg["logFC"] > 0.5)).sum()),
    down_strict=int(((deg["adjP"] < 0.05) & (deg["logFC"] < -0.5)).sum()),
)

# ---------------------------------------------------------------- 4. 18 genes
EPI = ["Krt12", "Aqp5", "Gsta4", "Upk3bl", "Upk1b", "Krt15", "Krt17",
       "Krt19", "Krt13", "Dsg1a", "Perp", "Tmprss11d", "Lypd3", "Calml3",
       "Gsta2", "Sfn", "Krt5", "Krt14"]
spec = pd.read_csv(r"第8天执行结果\results\D8_knockdown_spectrum.tsv", sep="\t",
                   usecols=["ko_gene", "affected_gene", "distance", "direction", "abs_distance"])
sel = spec[spec["ko_gene"].isin(["Tmsb4x", "Serpine1"]) & spec["affected_gene"].isin(EPI)]
tab = sel.pivot_table(index="affected_gene", columns="ko_gene", values="distance", aggfunc="first")
tab["abs_Tmsb4x"] = tab["Tmsb4x"].abs(); tab["abs_Serpine1"] = tab["Serpine1"].abs()
OUT["EPI18"] = {g: {c: (float(v) if pd.notna(v) else None) for c, v in row.items()}
                for g, row in tab.iterrows()}

json.dump(OUT, open("_rev12_stats.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)

# ---------------------------------------------------------------- print report
def f(x, n=4):
    return f"{x:.{n}g}" if isinstance(x, (int, float)) else str(x)

print("=== GSE241715 ===")
for g, d in stats241.items():
    print(f"{g:9s} Ctl {d['Control']['mean']:.3f}±{d['Control']['sd']:.3f} | KC {d['Keratoconus']['mean']:.3f} | "
          f"Vir {d['Viral']['mean']:.3f} | BF {d['Bacterial_fungal']['mean']:.3f} | "
          f"fold {f(d['fold_ctrl_to_bf'],3)} | AUC les {d['auc_lesion']:.3f} [{d['auc_lesion_ci'][0]:.3f},{d['auc_lesion_ci'][1]:.3f}] | "
          f"AUC BF {d['auc_bf']:.3f} [{d['auc_bf_ci'][0]:.3f},{d['auc_bf_ci'][1]:.3f}]")
print("LOOCV:", OUT["LOOCV"])
print("\n=== GSE58291 ===")
for g, d in res58291.items():
    print(f"{g:9s} {d['n']} means {d['mean']} | AUC BKvsCtl {d['auc_bk_vs_ctrl']:.3f} "
          f"[{d['auc_bk_vs_ctrl_ci'][0]:.3f},{d['auc_bk_vs_ctrl_ci'][1]:.3f}] | AUC InfvsCtl {d['auc_inf_vs_ctrl']:.3f}")
print("\n=== DEG ===", OUT["DEG"])
print("\n=== EPI18 (distance) ===")
print(tab.round(4).to_string())