# -*- coding: utf-8 -*-
"""图S4 多敲除叠加汇总
a  单敲与多敲的上皮 mean |Δ|：组合实测条形 + 单敲平均（空心菱形）+ 比值
b  实测 multi 对线性预测 individual_avg 的散点（y = x 虚线）

注：含 $...$ 的字符串整体交给 mathtext，中文字形会缺字，
   故含数学符号的文本一律只用 ASCII，中文文本一律不含 $。
"""
import numpy as np
import pandas as pd
import figstyle as S
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

df = S.load(9, "D9_multiko_summary.tsv")

ORDER = ["SERPINE1", "TMSB4X", "P4HB", "SERPINE1+TMSB4X", "SERPINE1+P4HB",
         "TMSB4X+P4HB", "SERPINE1+TMSB4X+P4HB"]
LAB = {"SERPINE1": "SERPINE1", "TMSB4X": "TMSB4X", "P4HB": "P4HB",
       "SERPINE1+TMSB4X": "S + T", "SERPINE1+P4HB": "S + P",
       "TMSB4X+P4HB": "T + P", "SERPINE1+TMSB4X+P4HB": "S + T + P（三敲）"}
df = df.set_index("combo").loc[ORDER].reset_index()

C1, C2, C3 = S.PAL["blue"], S.PAL["teal"], S.PAL["red"]
BAR_C = {1: C1, 2: C2, 3: C3}
print(df[["combo", "epi_mean_abs_multi", "epi_mean_abs_individual_avg",
          "ratio"]].to_string())

fig, axes = S.grid(1, 2, w=6.45, h=3.20, wspace=0.62, left=0.135, right=0.845,
                   top=0.885, bottom=0.255, widths=[1.0, 0.98])

# ---------------------------------------------------------------- (a)
ax = axes[0]
y = np.arange(len(ORDER))[::-1]
base = 1.15e-4
vals = df["epi_mean_abs_multi"].values
ind = df["epi_mean_abs_individual_avg"].values
cols = [BAR_C[int(n)] for n in df.n_genes]
ax.barh(y, vals - base, left=base, color=cols, height=0.60, edgecolor="none",
        zorder=2)
ax.scatter(ind, y, marker="D", s=15, facecolor="white", edgecolor="#333333",
           lw=0.7, zorder=4)
for yi, v, iv, r in zip(y, vals, ind, df["ratio"].values):
    ax.annotate(f"{v:.4f}", (v, yi), textcoords="offset points", xytext=(7, 0),
                va="center", ha="left", fontsize=5.6, color="#333333")
    ax.text(1.085, yi, f"×{r:.2f}", va="center", ha="left", fontsize=5.8,
            color="#333333", transform=ax.get_yaxis_transform(),
            clip_on=False, fontweight="bold" if r > 1.5 else "normal")
ax.set_xscale("log")
ax.set_xlim(base, 0.95)
ax.set_ylim(-0.62, len(ORDER) - 0.38)
ax.set_yticks(y)
ax.set_yticklabels([LAB[c] for c in ORDER], fontsize=6.2)
ax.set_xticks([1e-4, 1e-3, 1e-2, 1e-1])
ax.set_xticklabels(["$10^{-4}$", "$10^{-3}$", "$10^{-2}$", "$10^{-1}$"],
                   fontsize=6.0)
ax.set_xlabel("18 个上皮基因 mean |Δ|（对数轴）", fontsize=6.6, labelpad=3)
ax.set_title("组合敲除的实测效应", pad=6)
S.clean(ax)
ax.set_axisbelow(True)
ax.xaxis.grid(True, color="#EEEEEE", lw=0.5)
leg = [Patch(fc=C1, ec="none", label="单敲"),
       Patch(fc=C2, ec="none", label="双敲"),
       Patch(fc=C3, ec="none", label="三敲")]
l1 = ax.legend(handles=leg, loc="upper left", bbox_to_anchor=(0.0, -0.205),
               ncol=3, fontsize=5.4, handletextpad=0.45, borderpad=0.2,
               columnspacing=1.4, frameon=False)
ax.add_artist(l1)
leg2 = [Line2D([0], [0], marker="D", color="none", markerfacecolor="white",
               markeredgecolor="#333333", ms=4, label="单敲平均（线性预测）"),
        Line2D([0], [0], color="none", label="×ratio = 实测 / 单敲平均")]
ax.legend(handles=leg2, loc="upper left", bbox_to_anchor=(0.0, -0.29),
          ncol=2, fontsize=5.4, handletextpad=0.45, borderpad=0.2,
          columnspacing=1.4, frameon=False)
S.panel(ax, "a", x=-0.30)

# ---------------------------------------------------------------- (b)
ax = axes[1]
x = df["epi_mean_abs_individual_avg"].values
yv = df["epi_mean_abs_multi"].values
lo, hi = 3.0e-4, 0.85
gx = np.array([lo, hi])
ax.plot(gx, gx, ls=(0, (4, 2.5)), color="#8A8A8A", lw=0.8, zorder=1)
ax.scatter(x, yv, s=22, c=[BAR_C[int(n)] for n in df.n_genes],
           edgecolor="white", lw=0.6, zorder=4)

OFF = {
    "SERPINE1": (10, 5, "left"),
    "P4HB": (9, 3, "left"),
    "SERPINE1+P4HB": (-9, 6, "right"),
    "SERPINE1+TMSB4X": (0, 10, "center"),
    "TMSB4X+P4HB": (8, -9, "left"),
    "TMSB4X": (8, 2, "left"),
    "SERPINE1+TMSB4X+P4HB": (-8, -13, "right"),
}
NAME = {"SERPINE1": "SERPINE1", "P4HB": "P4HB", "SERPINE1+P4HB": "S + P",
        "SERPINE1+TMSB4X": "S + T", "TMSB4X+P4HB": "T + P",
        "TMSB4X": "TMSB4X", "SERPINE1+TMSB4X+P4HB": "三敲 ×2.81"}
for xi, yi, c in zip(x, yv, df.combo):
    o = OFF[c]
    ax.annotate(NAME[c], (xi, yi), textcoords="offset points",
                xytext=(o[0], o[1]), fontsize=5.6,
                color=S.PAL["red"] if c == "SERPINE1+TMSB4X+P4HB" else "#1A1A1A",
                ha=o[2], va="center")

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(lo, hi)
ax.set_ylim(lo, hi)
ax.set_aspect("equal")
tk = [1e-3, 1e-2, 1e-1]
ax.set_xticks(tk)
ax.set_yticks(tk)
ax.set_xticklabels(["$10^{-3}$", "$10^{-2}$", "$10^{-1}$"], fontsize=6.0)
ax.set_yticklabels(["$10^{-3}$", "$10^{-2}$", "$10^{-1}$"], fontsize=6.0)
ax.set_xlabel("线性预测：单敲平均 mean |Δ|", fontsize=6.6, labelpad=3)
ax.set_ylabel("实测组合 mean |Δ|", fontsize=6.6, labelpad=3)
ax.set_title("实测对线性预测", pad=6)
S.clean(ax)
ax.text(0.56, 0.055, "虚线：y = x（可加）", transform=ax.transAxes,
        fontsize=5.2, color="#777777", ha="left", va="bottom")
S.panel(ax, "b", x=-0.26)

S.save(fig, "S4_multiko")

