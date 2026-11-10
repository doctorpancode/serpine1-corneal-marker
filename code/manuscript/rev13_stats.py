# -*- coding: utf-8 -*-
"""v1.3 numbers: AUCell ordering, rank table, LOOCV for 10 genes, external calibration/DCA."""
import io
import json
import sys

import numpy as np
import pandas as pd

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
R = r"D:\A生信分析\虚拟敲除"
out = {}

# ------------------------------------------------------------------ 1. AUCell by celltype
auc = pd.read_csv(R + r"\第6天执行结果\results\D6_aucell_summary.tsv", sep="\t")
print("aucell genes:", sorted(auc.gene.unique()))
for g in ["Serpine1", "Tmsb4x"]:
    sub = auc[auc.gene == g]
    piv = sub.pivot_table(index="celltype", columns="sample", values="mean_auc")
    npiv = sub.pivot_table(index="celltype", columns="sample", values="n", aggfunc="sum")
    print(f"\n--- {g}: mean AUCell by celltype")
    t = piv.copy()
    t.columns = [f"auc_{c}" for c in t.columns]
    t["n_cells"] = npiv.sum(axis=1)
    print(t.sort_values([c for c in t.columns if c.startswith("auc_")][0], ascending=False).round(5).to_string())
    out[f"aucell_{g}"] = {k: {kk: (None if pd.isna(vv) else round(float(vv), 5)) for kk, vv in v.items()}
                          for k, v in t.to_dict("index").items()}

# ------------------------------------------------------------------ 2. rank table (4 metrics)
reg = pd.read_csv(R + r"\第9天执行结果\results\D9_regularized_ko_scores.tsv", sep="\t")
core = ["Tmsb4x", "Prdx2", "Capg", "Vdac1", "Xiap", "P4hb", "S100a13", "Sod1", "Serpine1", "Gja1", "Vim"]
tab = reg[reg.ko_gene.isin(core)][["ko_gene", "out_degree", "mean_raw", "mean_r1", "cos_to_tmsb4x",
                                   "rank_raw", "rank_r1", "rank_cos", "rank_resid"]].copy()
tab = tab.set_index("ko_gene").loc[core].reset_index()
print("\n--- rank table (4 metrics), n =", len(reg))
print(tab.to_string(index=False))
out["rank_table"] = tab.to_dict("records")
out["n_ko_total"] = int(len(reg))

# ------------------------------------------------------------------ 3. Krt12 rank distribution
spec = pd.read_csv(R + r"\第8天执行结果\results\D8_knockdown_spectrum.tsv", sep="\t")
k = spec[spec.affected_gene == "Krt12"].copy()
k["rank_in_ko"] = k.groupby("ko_gene").abs_distance.rank(ascending=False, method="min")
out["krt12_rank_median"] = float(k.rank_in_ko.median())
out["krt12_rank_q1"], out["krt12_rank_q3"] = [float(x) for x in k.rank_in_ko.quantile([0.25, 0.75])]
out["krt12_top10"] = int((k.rank_in_ko <= 10).sum())
print("\nKrt12 rank among genes perturbed by each KO: median", out["krt12_rank_median"],
      "IQR", (out["krt12_rank_q1"], out["krt12_rank_q3"]), "| top10 in", out["krt12_top10"], "KOs")
rows_per_ko = spec.groupby("ko_gene").size()
out["rows_per_ko_median"] = float(rows_per_ko.median())
print("rows per KO: median", out["rows_per_ko_median"], "min", int(rows_per_ko.min()), "max", int(rows_per_ko.max()))

# ------------------------------------------------------------------ 4. discovery AUC + LOOCV for 10 genes
expr = pd.read_csv(R + r"\D1执行结果\第一天执行结果\data\processed\GSE241715_expression.tsv", sep="\t")
sheet = pd.read_csv(R + r"\D1执行结果\第一天执行结果\data\processed\GSE241715_sample_sheet.tsv", sep="\t")
print("\nexpr:", expr.shape, "| sheet cols:", list(sheet.columns))
sym = expr.columns[0]
expr = expr.set_index(sym)
y = (sheet.set_index("sample").loc[expr.columns, "binary_lesion"].values).astype(int)
print("lesion:", int(y.sum()), "control:", int((1 - y).sum()))


def auc_score(y, p):
    from itertools import product
    pos = p[y == 1]
    neg = p[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    gt = sum((a > b) + 0.5 * (a == b) for a, b in product(pos, neg))
    return gt / (len(pos) * len(neg))


def logit_fit(x, y):
    X = np.column_stack([np.ones_like(x), x])
    b = np.zeros(2)
    for _ in range(200):
        p = 1 / (1 + np.exp(-(X @ b)))
        W = p * (1 - p) + 1e-10
        H = X.T @ (X * W[:, None])
        g = X.T @ (y - p)
        try:
            step = np.linalg.solve(H, g)
        except np.linalg.LinAlgError:
            break
        b = b + step
        if np.max(np.abs(step)) < 1e-10:
            break
    return b


genes = ["SERPINE1", "P4HB", "S100A13", "CAPG", "XIAP", "TMSB4X", "GJA1", "SOD1", "PRDX2", "VDAC1"]
res = []
for g in genes:
    if g not in expr.index:
        res.append({"gene": g, "auc": None, "loocv": None})
        continue
    x = expr.loc[g].values.astype(float)
    a = auc_score(y, x)
    # LOOCV: fit on n-1, predict held-out, pool
    pred = np.zeros(len(y))
    for i in range(len(y)):
        m = np.ones(len(y), bool)
        m[i] = False
        b = logit_fit(x[m], y[m])
        pred[i] = 1 / (1 + np.exp(-(b[0] + b[1] * x[i])))
    res.append({"gene": g, "auc": round(a, 4), "loocv": round(auc_score(y, pred), 4)})
print("\n--- discovery AUC / LOOCV")
print(pd.DataFrame(res).to_string(index=False))
out["discovery_auc_loocv"] = res

# ------------------------------------------------------------------ 5. external calibration + DCA
g58 = pd.read_csv(R + r"\第9天执行结果\data\GSE58291_expr.tsv", sep="\t")
grp = pd.read_csv(R + r"\第9天执行结果\results\D9_GSE58291_sample_group.tsv", sep="\t")
print("\nGSE58291 expr:", g58.shape, "| group cols:", list(grp.columns), grp.shape)
sym_col = g58.columns[1]
g58 = g58.set_index(sym_col)
import re
sm = {c: re.sub(r"\s+", "", c.split(".")[0]) for c in g58.columns if c.endswith("AVG_Signal")}
g58 = g58[[c for c in g58.columns if c.endswith("AVG_Signal")]]
g58.columns = [sm[c] for c in g58.columns]
print("samples:", len(g58.columns), "| group head:", grp.head(3).to_dict("records"))

gmap = dict(zip(grp.iloc[:, 0].astype(str), grp.iloc[:, 1].astype(str)))
keep = [c for c in g58.columns if c in gmap]
X = g58[keep]
lab = np.array([1 if gmap[c].upper().startswith(("BK", "FK", "B", "F")) and "CONTROL" not in gmap[c].upper() else 0 for c in keep])
print("kept:", len(keep), "| labels:", dict(zip(*np.unique(lab, return_counts=True))))

ser = X.loc["SERPINE1"].values.astype(float)
b = logit_fit(ser, lab)
p = 1 / (1 + np.exp(-(b[0] + b[1] * ser)))
lp = b[0] + b[1] * ser
# calibration slope/intercept: logistic on lp
bc = logit_fit(lp, lab)
out["calib_intercept"], out["calib_slope"] = float(bc[0]), float(bc[1])
out["ext_auc"] = round(auc_score(lab, ser), 4)
out["ext_auc_ci"] = None
# bootstrap CI
rng = np.random.default_rng(20261111)
bs = []
for _ in range(2000):
    idx = rng.integers(0, len(lab), len(lab))
    if len(np.unique(lab[idx])) < 2:
        continue
    bs.append(auc_score(lab[idx], ser[idx]))
out["ext_auc_ci"] = [round(float(np.percentile(bs, 2.5)), 3), round(float(np.percentile(bs, 97.5)), 3)]
# Hosmer-Lemeshow 3 groups (tertiles)
q = pd.qcut(p, 3, labels=False, duplicates="drop")
hl = 0.0
for g in np.unique(q):
    m = q == g
    o, e, n = lab[m].sum(), p[m].sum(), m.sum()
    if 0 < e < n:
        hl += (o - e) ** 2 / (e * (1 - e / n))
out["hl_chi2"], out["hl_df"] = float(hl), int(len(np.unique(q)) - 2)
from math import erf, sqrt
def chi2_sf(x, k):
    # survival function via regularized gamma (series for small x)
    if x <= 0:
        return 1.0
    a, s, term = k / 2, 0.0, 1.0 / (k / 2)
    # use numeric integration (simple, k<=3)
    n = 200000
    hi = 60.0
    h = hi / n
    tot = 0.0
    for i in range(1, n + 1):
        t = i * h
        tot += (t ** (k / 2 - 1)) * np.exp(-t / 2) * h
    from math import gamma
    return 1 - tot / (2 ** (k / 2) * gamma(k / 2))
out["hl_p"] = round(float(chi2_sf(hl, out["hl_df"])), 3)
# DCA: net benefit
prev = lab.mean()
dca = []
for th in np.arange(0.2, 0.81, 0.05):
    tp = ((p >= th) & (lab == 1)).sum()
    fp = ((p >= th) & (lab == 0)).sum()
    nb = tp / len(lab) - fp / len(lab) * th / (1 - th)
    nb_all = prev - (1 - prev) * th / (1 - th)
    dca.append({"threshold": round(float(th), 2), "nb_model": round(float(nb), 3),
                "nb_all": round(float(nb_all), 3), "nb_none": 0.0})
print("\n--- external calibration: AUC", out["ext_auc"], out["ext_auc_ci"],
      "| slope", round(out["calib_slope"], 3), "intercept", round(out["calib_intercept"], 3),
      "| HL chi2", round(out["hl_chi2"], 2), "df", out["hl_df"], "P", out["hl_p"])
print(pd.DataFrame(dca).to_string(index=False))
out["dca"] = dca
out["ext_prev"] = round(float(prev), 3)

# ------------------------------------------------------------------ 6. immune counts
imm = pd.read_csv(R + r"\第5天执行结果\results\D5_immune_group_stats.tsv", sep="\t")
out["immune_n_sig"] = int((imm.adjP < 0.05).sum())
print("\nimmune celltypes adjP<0.05:", out["immune_n_sig"], "of", len(imm))
print(imm.sort_values("adjP").head(6)[["celltype", "delta", "p_wilcox", "adjP"]].to_string(index=False))

json.dump(out, open(R + r"\_rev13_stats.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1, default=str)
print("\nwritten _rev13_stats.json")