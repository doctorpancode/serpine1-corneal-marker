# -*- coding: utf-8 -*-
"""Figure S10 — rank stability across 20 network rebuilds and 20 degree-preserving nulls."""
import io
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

reb = pd.read_csv(os.path.join(RES, "D9_rank_stability.tsv"), sep="\t")
nul = pd.read_csv(os.path.join(RES, "D9_null_ranks.tsv"), sep="\t")

fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.5),
                         gridspec_kw={"width_ratios": [2.05, 1]})

# ---- (a) strip plot: rebuild vs null, residual ranks 1..947 (log scale)
ax = axes[0]
rng = np.random.default_rng(7)
ys = np.arange(len(CORE))[::-1]
for j, g in enumerate(CORE):
    for data, off, col, lab in ((reb, 0.17, "#B2182B", "重建"), (nul, -0.17, "#4A6FA5", "零分布")):
        v = data.loc[data.ko_gene == g, "rank_resid"].values
        x = np.clip(v, 1, 947)
        ax.scatter(x, np.full_like(x, ys[j] + off, dtype=float) + rng.uniform(-0.05, 0.05, len(x)),
                   s=5.5, color=col, alpha=0.55, linewidths=0,
                   label=lab if j == 0 else None)
ax.set_xscale("log")
ax.set_yticks(ys)
ax.set_yticklabels([LAB[g] for g in CORE], fontsize=6.8)
ax.set_xlabel("上皮残差排名（1 = 最高，947 = 最低，log 轴）")
ax.set_xlim(0.85, 1100)
ax.axvline(947 / 3, color="0.55", lw=0.6, ls=(0, (3, 2)))
ax.text(947 / 3 * 1.12, len(CORE) - 0.45, "后 2/3 起点", fontsize=6, color="0.35")
ax.legend(loc="lower left", fontsize=6.2, frameon=False, ncol=1,
          handletextpad=0.3, columnspacing=0.8)
ax.set_title("a  20 次重建与 20 次度保持置换", fontsize=7.6, loc="left", pad=4)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

# ---- (b) bar: share of draws in the bottom third
ax = axes[1]
share = []
for g in CORE:
    r = (reb.loc[reb.ko_gene == g, "rank_resid"] > 947 / 3).mean() * 100
    n = (nul.loc[nul.ko_gene == g, "rank_resid"] > 947 / 3).mean() * 100
    share.append((r, n))
share = np.array(share)
y = np.arange(len(CORE))[::-1]
ax.barh(y + 0.18, share[:, 0], height=0.34, color="#B2182B", label="重建")
ax.barh(y - 0.18, share[:, 1], height=0.34, color="#4A6FA5", label="零分布")
ax.set_yticks(y)
ax.set_yticklabels([])
ax.set_xlim(0, 105)
ax.set_xlabel("落于后三分之一的次数占比（%）")
ax.legend(loc="upper left", bbox_to_anchor=(0.02, 0.99), fontsize=6.2, frameon=False)
ax.set_title("b  位置稳定性", fontsize=7.6, loc="left", pad=4)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

fig.tight_layout(pad=0.4)
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(OUT, f"FigS10_rank_stability.{ext}"), dpi=400,
                bbox_inches="tight")
print("written", os.path.join(OUT, "FigS10_rank_stability.png"))