# -*- coding: utf-8 -*-
"""图S6 —— 多敲除逐基因效应矩阵（18 个上皮基因 × 3 个单敲）"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

import figstyle as S

d = S.load(9, "D9_multiko_epi_table.tsv")

COLS = ["SERPINE1", "TMSB4X", "P4HB"]
mat = d[COLS].to_numpy(float)
# 列内归一化（各列量级相差 3 个数量级，共用色标会把 SERPINE1 压成白色）
scale = np.abs(mat).max(axis=0)
z = mat / np.where(scale == 0, 1, scale)
abs_sum = d["abs_sum"].to_numpy(float)
abs3 = np.abs(mat).sum(axis=1)          # 三列绝对值之和（含 P4HB）
same = d["same_sign"].to_numpy(bool)
genes = list(d["gene"])
n = len(genes)

fig, axes = S.grid(1, 2, w=6.3, h=4.75, widths=[1.05, 1.0],
                   wspace=0.50, left=0.135, right=0.985, top=0.905,
                   bottom=0.255)

# ---------------------------------------------------------------- (a)
ax = axes[0]
annot = [["" for _ in COLS] for _ in range(n)]
for i in range(n):
    for j in range(len(COLS)):
        v = mat[i, j]
        annot[i][j] = f"{v:.2f}" if abs(v) >= 1e-2 else f"{v:.1e}"
S.heat(ax, z, genes, COLS, cmap=S.CMAP_DIV, vmin=-1, vmax=1, center=0,
       cbar_label="Δ 表达（列内归一化）", annot=annot, annot_size=4.9,
       xrot=30, cbar_kw=dict(fraction=0.042, pad=0.03))
ax.set_title("三敲逐基因 Δ 表达", pad=6)
ax.tick_params(axis="y", labelsize=6.0)
ax.tick_params(axis="x", labelsize=6.0)
S.panel(ax, "a", x=-0.46, y=1.045)

# ---------------------------------------------------------------- (b)
ax = axes[1]
y = np.arange(n)
col = [S.PAL["teal"] if s else S.PAL["brick"] for s in same]
ax.barh(y, abs_sum, height=0.66, color=col, edgecolor="white", lw=0.3,
        zorder=3)
ax.scatter(abs3, y, marker="|", s=26, color="#3A3A3A", lw=0.8, zorder=4)

for yi, v, v3 in zip(y, abs_sum, abs3):
    ax.annotate(f"{v:.3f}", (max(v, v3), yi), textcoords="offset points", xytext=(4, 0),
                va="center", ha="left", fontsize=4.9, color="#333333")

ax.set_ylim(n - 0.45, -0.55)
ax.set_xlim(0, 0.90)
ax.set_yticks(y)
ax.set_yticklabels(genes, fontsize=5.8)
ax.tick_params(axis="y", length=2, width=0.5, pad=1.5)
ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8])
ax.tick_params(axis="x", labelsize=6.0)
ax.set_xlabel("三敲绝对效应 abs_sum", fontsize=6.6, labelpad=3)
ax.set_title("效应量与符号一致性", pad=6)
S.clean(ax)
for yi in y:
    if yi % 2 == 1:
        ax.axhspan(yi - 0.5, yi + 0.5, color="#F7F7F7", lw=0, zorder=0)

leg = [Patch(fc=S.PAL["teal"], ec="white",
             label=f"ΔSERPINE1 与 ΔTMSB4X 同号（{int(same.sum())}/{n}）"),
       Patch(fc=S.PAL["brick"], ec="white",
             label=f"两者异号（{int((~same).sum())}/{n}）"),
       plt.Line2D([0], [0], marker="|", ls="none", color="#3A3A3A",
                  ms=6, mew=0.9, label="三列绝对值之和")]
ax.legend(handles=leg, loc="upper center", bbox_to_anchor=(0.5, -0.085),
          ncol=1, fontsize=5.0, frameon=False, handletextpad=0.5,
          handlelength=1.1, labelspacing=0.3, borderpad=0.1)
ax.get_legend().remove()
fig.legend(handles=leg, loc="lower center", bbox_to_anchor=(0.5, 0.128),
           ncol=3, fontsize=5.2, frameon=False, handletextpad=0.5,
           handlelength=1.1, columnspacing=1.6, borderpad=0.1)
S.panel(ax, "b", x=-0.26, y=1.045)

fig.text(0.135, 0.075,
         "行序按 |ΔTMSB4X| 降序；热图色标 = 各列值 / 该列最大绝对值，"
         "格内数字为原始 Δ 表达（对数记法）",
         fontsize=5.2, color="#666666", ha="left", va="bottom")
fig.text(0.135, 0.030,
         "abs_sum = |ΔSERPINE1| + |ΔTMSB4X|（表内定义，未含 P4HB）；"
         "细线 = 三列绝对值之和；同号判定依据 ΔSERPINE1 与 ΔTMSB4X 的符号",
         fontsize=5.2, color="#666666", ha="left", va="bottom")

S.save(fig, "S6_multiko_genes")
plt.close(fig)

