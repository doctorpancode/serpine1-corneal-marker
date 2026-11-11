# -*- coding: utf-8 -*-
"""Figure S11 — effect magnitude tracks expression abundance; class comparison."""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = r"D:\A生信分析\虚拟敲除"
RES = os.path.join(BASE, "第9天执行结果", "results")
OUT = os.path.join(BASE, "文章", "figs_v2", "out")
CORE = ["Tmsb4x", "Prdx2", "Capg", "Vdac1", "P4hb", "Xiap", "S100a13", "Sod1",
        "Serpine1", "Vim", "Gja1"]
LAB = {"Tmsb4x": "TMSB4X", "Prdx2": "PRDX2", "Capg": "CAPG", "Vdac1": "VDAC1",
       "P4hb": "P4HB", "Xiap": "XIAP", "S100a13": "S100A13", "Sod1": "SOD1",
       "Serpine1": "SERPINE1", "Vim": "VIM", "Gja1": "GJA1"}

plt.rcParams.update({"font.family": ["Arial", "Microsoft YaHei"], "font.size": 7.2,
                     "axes.unicode_minus": False,
                     "axes.linewidth": 0.7, "xtick.major.width": 0.6,
                     "ytick.major.width": 0.6, "pdf.fonttype": 42})

d = pd.read_csv(os.path.join(RES, "D9_control_analyses.tsv"), sep="\t", index_col=0)
# fraction of the same-decile genes lying below, in percent (consistent with the text)
d["pct_below"] = [100.0 * (d.loc[d["decile"] == d.loc[g, "decile"], "mean_raw"]
                           < d.loc[g, "mean_raw"]).mean() for g in d.index]

fig, axes = plt.subplots(1, 3, figsize=(6.9, 2.45),
                         gridspec_kw={"width_ratios": [1.35, 1.0, 0.95]})

# ---- (a) effect vs expression
ax = axes[0]
floor = 3e-5
y = d["mean_raw"].clip(lower=floor)
ax.scatter(d["expr"], y, s=3.2, color="0.74", linewidths=0, alpha=0.7)
for g in CORE:
    ax.scatter(d.loc[g, "expr"], max(d.loc[g, "mean_raw"], floor), s=12,
               color="#B2182B" if g == "Tmsb4x" else "#4A6FA5", linewidths=0, zorder=3)
top = d.nlargest(3, "mean_raw").index.tolist()
ax.scatter(d.loc[top, "expr"], d.loc[top, "mean_raw"], s=22, facecolors="none",
           edgecolors="#B2182B", linewidths=0.8, zorder=4)
ax.annotate("TMSB4X", (d.loc["Tmsb4x", "expr"], d.loc["Tmsb4x", "mean_raw"]),
            xytext=(-30, 8), textcoords="offset points", fontsize=6.4, color="0.12",
            arrowprops=dict(arrowstyle="-", lw=0.5, color="0.45"))
ax.annotate("mt-Co1", (d.loc["mt-Co1", "expr"], d.loc["mt-Co1", "mean_raw"]),
            xytext=(6, -13), textcoords="offset points", fontsize=6.4, color="0.12")
ax.annotate("S100A11", (d.loc["S100a11", "expr"], d.loc["S100a11", "mean_raw"]),
            xytext=(-14, -16), textcoords="offset points", fontsize=6.4, color="0.12")
ax.annotate("SERPINE1", (d.loc["Serpine1", "expr"], d.loc["Serpine1", "mean_raw"]),
            xytext=(10, 16), textcoords="offset points", fontsize=6.4, color="0.12",
            arrowprops=dict(arrowstyle="-", lw=0.5, color="0.45"))
ax.set_yscale("log")
ax.set_ylim(floor * 0.75, 0.85)
ax.set_xlim(-0.3, 11.4)
ax.set_xlabel("基因平均表达（log2 反卷积签名，11 类细胞）")
ax.set_ylabel("上皮平均核效应 mean|Δ|")
ax.set_title("a  效应量随表达丰度上升", fontsize=7.6, loc="left", pad=4)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

# ---- (b) within-decile percentile for the 11 core genes
ax = axes[1]
ys = np.arange(len(CORE))[::-1]
pct = [d.loc[g, "pct_below"] for g in CORE]
cols = ["#B2182B" if p < 50 else "#4A6FA5" for p in pct]
ax.barh(ys, pct, height=0.62, color=cols)
ax.axvline(50, color="0.4", lw=0.6, ls=(0, (3, 2)))
ax.set_yticks(ys)
ax.set_yticklabels([LAB[g] for g in CORE], fontsize=6.8)
ax.set_xlim(0, 118)
ax.set_xlabel("同表达十分位内的百分位（%，50 = 层内中位）")
ax.set_title("b  丰度匹配后的位置", fontsize=7.6, loc="left", pad=4)
for yy, p in zip(ys, pct):
    ax.text(p + 2.5, yy, f"{p:.0f}", va="center", fontsize=6.0, color="0.25")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

# ---- (c) secreted vs transcription factor
ax = axes[2]
order = [("secreted", "分泌型"), ("tf", "转录因子"), ("neither", "其余")]
rng = np.random.default_rng(11)
for i, (k, lab) in enumerate(order):
    v = d.loc[d["class"] == k, "mean_raw"].clip(lower=floor)
    ax.scatter(np.full(len(v), i) + rng.uniform(-0.13, 0.13, len(v)), v, s=2.6,
               color="0.62", linewidths=0, alpha=0.65)
    ax.plot([i - 0.28, i + 0.28], [v.median()] * 2, color="#B2182B", lw=1.1, zorder=3)
    ax.text(i, 0.30, f"n={len(v)}", ha="center", fontsize=6.0, color="0.3")
ax.set_yscale("log")
ax.set_ylim(floor * 0.75, 0.60)
ax.set_xticks(range(3))
ax.set_xticklabels([lab for _, lab in order], fontsize=6.8)
ax.set_yticks([1e-4, 1e-3, 1e-2, 1e-1])
ax.set_ylabel("mean|Δ|")
ax.set_title("c  类别间无差异", fontsize=7.6, loc="left", pad=4)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

fig.tight_layout(pad=0.4)
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(OUT, f"FigS11_expression_confound.{ext}"), dpi=400,
                bbox_inches="tight")
print("written", os.path.join(OUT, "FigS11_expression_confound.png"))