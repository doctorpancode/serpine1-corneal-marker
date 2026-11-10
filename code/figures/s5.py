# -*- coding: utf-8 -*-
"""图S5 归一化前后的影响谱重排
a  947 个敲除按 mean_raw 降序的 mean |Δ| 谱（对数纵轴），11 个核心基因高亮
b  11 个核心基因的排名变化 dumbbell（rank_raw -> rank_r1）

注：含 $...$ 的字符串整体交给 mathtext，中文字形会缺字，
   故含数学符号的文本一律只用 ASCII，中文文本一律不含 $。
"""
import numpy as np
import pandas as pd
import figstyle as S
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

r = S.load(9, "D9_regularized_ko_scores.tsv")
r = r.sort_values("rank_raw").reset_index(drop=True)

CORE = ["Tmsb4x", "Prdx2", "Capg", "P4hb", "Gja1", "Vdac1", "S100a13",
        "Vim", "Sod1", "Serpine1", "Xiap"]
print(r[r.ko_gene.isin(CORE)][["ko_gene", "mean_raw", "rank_raw", "rank_r1"]]
      .to_string())

fig, axes = S.grid(1, 2, w=6.55, h=3.45, wspace=0.50, left=0.085, right=0.975,
                   top=0.895, bottom=0.155, widths=[1.18, 1.0])

# ---------------------------------------------------------------- (a)
ax = axes[0]
x = np.arange(1, len(r) + 1)
y0 = 1.0e-7
zo = r["mean_raw"] < y0
ax.scatter(x[~zo], r.loc[~zo, "mean_raw"], s=1.3, c="#C9C9C9", lw=0,
           rasterized=True, zorder=2)
ax.scatter(x[zo], np.full(zo.sum(), y0 * 1.25), s=1.6, c="#B0B0B0",
           marker="v", lw=0, rasterized=True, zorder=2)

LABOFF = {
    "Tmsb4x": (8, 4, "left"), "Prdx2": (8, 7, "left"),
    "Capg": (10, -34, "left"), "P4hb": (-8, -11, "right"),
    "Gja1": (8, 9, "left"), "Vdac1": (8, -10, "left"),
    "S100a13": (8, 7, "left"), "Vim": (8, -11, "left"),
    "Sod1": (8, 9, "left"), "Serpine1": (8, 7, "left"),
    "Xiap": (8, -8, "left"),
}
for g in CORE:
    row = r[r.ko_gene == g].iloc[0]
    ax.scatter([row.rank_raw], [row.mean_raw], s=17, facecolor=S.PAL["red"],
               edgecolor="white", lw=0.7, zorder=5)
    dx, dy, ha = LABOFF[g]
    if g in ("Capg", "P4hb"):     # 左侧两点标签远离数据点，加引线
        ax.annotate("", xy=(row.rank_raw, row.mean_raw), xytext=(dx, dy),
                    textcoords="offset points", zorder=4,
                    arrowprops=dict(arrowstyle="-", lw=0.45, color="#9A9A9A",
                                    shrinkA=0.5, shrinkB=2.0))
    ax.annotate(f"{g}", (row.rank_raw, row.mean_raw), textcoords="offset points",
                xytext=(dx, dy), fontsize=5.3, color="#1A1A1A", ha=ha,
                va="center", zorder=6,
                bbox=dict(fc="white", ec="none", pad=0.7, alpha=0.85))

ax.set_yscale("log")
ax.set_xlim(0, 960)
ax.set_ylim(y0, 1.4)
ax.set_yticks([1e-7, 1e-5, 1e-3, 1e-1])
ax.set_yticklabels(["$10^{-7}$", "$10^{-5}$", "$10^{-3}$", "$10^{-1}$"],
                   fontsize=5.8)
ax.set_xticks([1, 200, 400, 600, 800, 947])
ax.set_xticklabels(["1", "200", "400", "600", "800", "947"], fontsize=6.0)
ax.set_xlabel("敲除按 mean |Δ| 降序排列的名次", fontsize=6.6, labelpad=3)
ax.set_ylabel("mean |Δ|（对数轴）", fontsize=6.6, labelpad=3)
ax.set_title("947 个敲除的影响谱", pad=6)
S.clean(ax)
leg = [Line2D([0], [0], marker="o", color="none", markerfacecolor=S.PAL["red"],
              markeredgecolor="white", ms=3.6, label="11 个核心基因"),
       Line2D([0], [0], marker="o", color="none", markerfacecolor="#C9C9C9",
              markeredgecolor="none", ms=2.2, label="其他敲除")]
ax.legend(handles=leg, loc="upper right", fontsize=5.3, handletextpad=0.4,
          borderpad=0.25, labelspacing=0.28)
ax.text(0.985, 0.30, "排名越靠后，影响越弱", transform=ax.transAxes,
        ha="right", va="top", fontsize=5.2, color="#888888")
ax.text(0.02, 0.015, f"下缘三角：{int(zo.sum())} 个敲除 mean |Δ| < 1e-7",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=5.0,
        color="#999999")
S.panel(ax, "a", x=-0.13)

# ---------------------------------------------------------------- (b)
ax = axes[1]
d = r[r.ko_gene.isin(CORE)].sort_values("rank_raw").reset_index(drop=True)
y = np.arange(len(d))[::-1]
for yi, rr in zip(y, d.itertuples()):
    up = rr.rank_r1 < rr.rank_raw
    c = S.PAL["blue"] if up else (S.PAL["red"] if rr.rank_r1 > rr.rank_raw
                                  else "#8C8C8C")
    ax.plot([rr.rank_raw, rr.rank_r1], [yi, yi], color=c, lw=1.0, zorder=2,
            solid_capstyle="round")
    ax.scatter([rr.rank_raw], [yi], s=15, facecolor="white", edgecolor=c,
               lw=0.8, zorder=3)
    ax.scatter([rr.rank_r1], [yi], s=17, facecolor=c, edgecolor="white",
               lw=0.6, zorder=4)
    hi, lo = max(rr.rank_raw, rr.rank_r1), min(rr.rank_raw, rr.rank_r1)
    ax.annotate(f"{rr.rank_raw}→{rr.rank_r1}", (hi, yi), textcoords="offset points",
                xytext=(6, 0), va="center", ha="left", fontsize=5.4,
                color="#333333")

ax.set_yticks(y)
ax.set_yticklabels(d.ko_gene, fontsize=6.4)
ax.set_xlim(0, 1000)         # 排名 1 在左侧（名次最优侧）
ax.set_ylim(-0.65, len(d) - 0.35)
ax.set_xticks([1, 200, 400, 600, 800, 947])
ax.set_xticklabels(["1", "200", "400", "600", "800", "947"], fontsize=6.0)
ax.set_xlabel("网络影响排名（1 = 最强，位于左侧）", fontsize=6.6, labelpad=3)
ax.set_title("行 L1 归一化前后的排名位移", pad=6)
S.clean(ax)
ax.set_axisbelow(True)
ax.xaxis.grid(True, color="#EEEEEE", lw=0.5)
leg = [Line2D([0], [0], marker="o", color="none", markerfacecolor="white",
              markeredgecolor="#555555", ms=3.6, label="归一化前 rank_raw"),
       Line2D([0], [0], marker="o", color="none", markerfacecolor="#555555",
              markeredgecolor="white", ms=3.6, label="归一化后 rank_r1"),
       Line2D([0], [0], color=S.PAL["blue"], lw=1.0, label="排名上升"),
       Line2D([0], [0], color=S.PAL["red"], lw=1.0, label="排名下降")]
ax.legend(handles=leg, loc="lower left", bbox_to_anchor=(-0.02, -0.30),
          ncol=2, fontsize=5.3, handletextpad=0.45, columnspacing=1.0,
          borderpad=0.2, labelspacing=0.3, frameon=False)
S.panel(ax, "b", x=-0.30)

S.save(fig, "S5_spectrum_shift")


