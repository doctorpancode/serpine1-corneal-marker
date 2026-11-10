# -*- coding: utf-8 -*-
"""D9_GSE58291_TMSB4X.tsv — both TMSB4X probes in GSE58291, per-sample, with AUC per probe."""
import io
import re

import numpy as np
import pandas as pd

EXPR = r"D:\A生信分析\虚拟敲除\第9天执行结果\data\GSE58291_expr.tsv"
GRP = r"D:\A生信分析\虚拟敲除\第9天执行结果\results\D9_GSE58291_sample_group.tsv"
OUT = r"D:\A生信分析\虚拟敲除\第9天执行结果\results\D9_GSE58291_TMSB4X.tsv"

expr = pd.read_csv(EXPR, sep="\t", low_memory=False)
cols = [c for c in expr.columns if isinstance(c, str) and c.endswith(".AVG_Signal")]
mapper = {c: re.sub(r"\s+", "", c.split(".")[0]) for c in cols}
grp = pd.read_csv(GRP, sep="\t")
gmap = grp.set_index("sample")["group"]


def auc(y, x):
    pos, neg = x[y == 1], x[y == 0]
    return float(np.mean([(p > q) + 0.5 * (p == q) for p in pos for q in neg]))


def boot_ci(y, x, n=2000, seed=20261111):
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        idx = rng.integers(0, len(y), len(y))
        if len(set(y[idx])) < 2:
            continue
        out.append(auc(y[idx], x[idx]))
    return float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


rows = []
for _, r in expr[expr.SYMBOL == "TMSB4X"].iterrows():
    s = pd.Series({mapper[c]: r[c] for c in cols})
    s = s[s.index.isin(gmap.index)]
    for samp, val in s.items():
        rows.append({"probe": r["PROBE_ID"], "sample": samp, "group": gmap[samp],
                     "signal": float(val), "log2_signal": float(np.log2(max(val, 1)))})
out = pd.DataFrame(rows).sort_values(["probe", "group", "sample"])
out.to_csv(OUT, sep="\t", index=False, float_format="%.4f")

for probe, d in out.groupby("probe"):
    y = (d["group"] != "Control").astype(int).values
    x = d["signal"].values
    a = auc(y, x)
    lo, hi = boot_ci(y, x)
    m = d.groupby("group")["signal"].mean().round(1).to_dict()
    print(f"{probe}  n={len(d)}  means={m}  AUC={a:.4f} ({lo:.3f}-{hi:.3f})")
print("rows", len(out), "->", OUT)