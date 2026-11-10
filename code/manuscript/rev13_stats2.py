# -*- coding: utf-8 -*-
"""v1.3 follow-up: spectrum sparsity/ties, transported calibration + DCA, immune counts."""
import io
import json
import re
import sys

import numpy as np
import pandas as pd

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
R = r"D:\A生信分析\虚拟敲除"
out = {}

# -------------------------------------------------- spectrum: how sparse / how tied
spec = pd.read_csv(R + r"\第8天执行结果\results\D8_knockdown_spectrum.tsv", sep="\t")
n = len(spec)
z = int((spec.abs_distance == 0).sum())
tiny = int((spec.abs_distance < 1e-12).sum())
out["spectrum_rows"] = n
out["frac_zero"] = round(z / n, 4)
out["frac_lt1e-12"] = round(tiny / n, 4)
print(f"spectrum rows {n}: exact zeros {z} ({z/n:.1%}), <1e-12 {tiny}")

# per-KO number of distinct non-zero affected genes
g = spec[spec.abs_distance > 1e-12].groupby("ko_gene").size()
out["nonzero_per_ko_median"] = float(g.median())
out["nonzero_per_ko_min"] = int(g.min())
out["nonzero_per_ko_max"] = int(g.max())
print("non-zero affected genes per KO: median", g.median(), "range", g.min(), "-", g.max())

# Krt12: rank among non-zero genes only, and its own |delta|
k = spec[spec.affected_gene == "Krt12"][["ko_gene", "abs_distance", "direction"]]
k = k.merge(g.rename("n_nonzero"), left_on="ko_gene", right_index=True)
nz = spec[spec.abs_distance > 1e-12].copy()
nz["rk"] = nz.groupby("ko_gene").abs_distance.rank(ascending=False, method="first")
kr = nz[nz.affected_gene == "Krt12"][["ko_gene", "rk", "abs_distance"]]
out["krt12_nonzero_median_rank"] = float(kr.rk.median())
out["krt12_nonzero_top10"] = int((kr.rk <= 10).sum())
out["krt12_zero_ko"] = int((k.abs_distance <= 1e-12).sum())
print("Krt12: median rank among non-zero genes", kr.rk.median(),
      "| top10 in", int((kr.rk <= 10).sum()), "KOs | zero delta in", out["krt12_zero_ko"], "KOs")
print("Krt12 |delta| quantiles:", [round(float(x), 4) for x in k.abs_distance.quantile([.25, .5, .75, .95])])

# -------------------------------------------------- external: transported calibration
expr = pd.read_csv(R + r"\D1执行结果\第一天执行结果\data\processed\GSE241715_expression.tsv", sep="\t").set_index("SYMBOL")
sheet = pd.read_csv(R + r"\D1执行结果\第一天执行结果\data\processed\GSE241715_sample_sheet.tsv", sep="\t")
y1 = sheet.set_index("sample").loc[expr.columns, "binary_lesion"].values.astype(int)
x1 = expr.loc["SERPINE1"].values.astype(float)


def logit_fit(X, y, iters=300):
    b = np.zeros(X.shape[1])
    for _ in range(iters):
        eta = np.clip(X @ b, -30, 30)
        p = 1 / (1 + np.exp(-eta))
        W = np.clip(p * (1 - p), 1e-9, None)
        H = X.T @ (X * W[:, None])
        gr = X.T @ (y - p)
        try:
            step = np.linalg.solve(H, gr)
        except np.linalg.LinAlgError:
            break
        b = b + step
        if np.max(np.abs(step)) < 1e-9:
            break
    return b


def auc(y, s):
    pos, neg = s[y == 1], s[y == 0]
    return float(np.mean([(a > b) + 0.5 * (a == b) for a in pos for b in neg]))


z1 = (x1 - x1.mean()) / x1.std(ddof=1)
X1 = np.column_stack([np.ones_like(z1), z1])
b1 = logit_fit(X1, y1)
out["disc_b0"], out["disc_b1"] = float(b1[0]), float(b1[1])
print("\ndiscovery logit on z(SERPINE1): b0 =", round(b1[0], 4), " b1 =", round(b1[1], 4))

g58 = pd.read_csv(R + r"\第9天执行结果\data\GSE58291_expr.tsv", sep="\t", low_memory=False).set_index("SYMBOL")
grp = pd.read_csv(R + r"\第9天执行结果\results\D9_GSE58291_sample_group.tsv", sep="\t")
g58 = g58[[c for c in g58.columns if str(c).endswith("AVG_Signal")]]
g58.columns = [re.sub(r"\s+", "", str(c).split(".")[0]) for c in g58.columns]
gmap = dict(zip(grp.iloc[:, 0].astype(str), grp.iloc[:, 1].astype(str)))
keep = [c for c in g58.columns if c in gmap]
x2 = g58.loc["SERPINE1", keep].values.astype(float)
y2 = np.array([0 if gmap[c].lower().startswith("cont") else 1 for c in keep])
print("external n =", len(y2), "| infected", int(y2.sum()), "| control", int((1 - y2).sum()))

z2 = (x2 - x2.mean()) / x2.std(ddof=1)
lp2 = b1[0] + b1[1] * z2
p2 = 1 / (1 + np.exp(-lp2))
out["ext_auc_transported"] = round(auc(y2, lp2), 4)
rng = np.random.default_rng(20261111)
bs = [auc(y2[i], lp2[i]) for i in (rng.integers(0, len(y2), len(y2)) for _ in range(2000))
      if len(np.unique(y2[i])) == 2]
out["ext_auc_ci"] = [round(float(np.percentile(bs, 2.5)), 3), round(float(np.percentile(bs, 97.5)), 3)]
print("transported AUC", out["ext_auc_transported"], out["ext_auc_ci"])

# calibration slope / intercept of the transported model (Cox calibration)
bc = logit_fit(np.column_stack([np.ones_like(lp2), lp2]), y2)
out["calib_slope"], out["calib_intercept"] = round(float(bc[1]), 3), round(float(bc[0]), 3)
print("calibration slope", out["calib_slope"], "intercept", out["calib_intercept"])

# calibration-in-the-large: observed vs mean predicted
out["obs_rate"], out["mean_pred"] = round(float(y2.mean()), 3), round(float(p2.mean()), 3)
print("observed", out["obs_rate"], "mean predicted", out["mean_pred"])

# HL, 3 groups
qq = pd.qcut(p2, 3, labels=False, duplicates="drop")
hl = 0.0
for gg in np.unique(qq):
    m = qq == gg
    o, e, nn = y2[m].sum(), p2[m].sum(), m.sum()
    hl += (o - e) ** 2 / (e * (1 - e / nn))
from math import gamma
def chi2_sf(x, k):
    nsteps = 200000
    hi, h = 60.0, 60.0 / 200000
    tot = sum(((i * h) ** (k / 2 - 1)) * np.exp(-i * h / 2) * h for i in range(1, nsteps + 1))
    return float(1 - tot / (2 ** (k / 2) * gamma(k / 2)))
out["hl_chi2"], out["hl_df"] = round(float(hl), 3), int(len(np.unique(qq)) - 2)
out["hl_p"] = round(chi2_sf(hl, out["hl_df"]), 3)
print("HL chi2", out["hl_chi2"], "df", out["hl_df"], "P", out["hl_p"])

# DCA
prev = float(y2.mean())
rows = []
for th in np.arange(0.1, 0.91, 0.05):
    flag = p2 >= th
    tp, fp = int((flag & (y2 == 1)).sum()), int((flag & (y2 == 0)).sum())
    nb = tp / len(y2) - fp / len(y2) * th / (1 - th)
    nb_all = prev - (1 - prev) * th / (1 - th)
    rows.append({"threshold": round(float(th), 2), "nb_model": round(float(nb), 3), "nb_all": round(float(nb_all), 3)})
dca = pd.DataFrame(rows)
print(dca.to_string(index=False))
out["dca"] = rows
out["nb_range"] = [round(float(dca.nb_model.max()), 3), round(float(dca.loc[dca.nb_model > dca.nb_all, "threshold"].min()), 2)]

# -------------------------------------------------- immune
imm = pd.read_csv(R + r"\第5天执行结果\results\D5_immune_group_stats.tsv", sep="\t")
out["immune_n_sig"] = int((imm.adjP < 0.05).sum())
print("\nimmune adjP<0.05:", out["immune_n_sig"], "of", len(imm))
print(imm.sort_values("adjP").head(6)[["celltype", "delta", "adjP"]].to_string(index=False))

json.dump(out, open(R + r"\_rev13_stats2.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1, default=str)
print("\nwritten _rev13_stats2.json")

