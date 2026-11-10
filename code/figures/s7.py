# -*- coding: utf-8 -*-
"""图S7  GSE326016 上皮分化程序得分（DM vs 正常，n = 3 / 3）
a 两组散点 + 均值/SD；b 按排名配对的 dumbbell（逐样本看方向）
统计为 Welch t 检验（对应稿件 P = 0.1097；等方差 t 检验 P = 0.0857）。
"""
import numpy as np
import pandas as pd
from scipy import stats
import figstyle as S

e = S.load(9, "D9_GSE326016_epi_score.tsv")
dm = e.loc[e.group == "DM", "epi_score"].to_numpy(float)
nm = e.loc[e.group == "Normal", "epi_score"].to_numpy(float)
t_w, p_w = stats.ttest_ind(dm, nm, equal_var=False)
print("DM", dm.mean(), dm.std(ddof=1), "NM", nm.mean(), nm.std(ddof=1))
print("welch t=%.3f p=%.4f" % (t_w, p_w),
      "student p=%.4f" % stats.ttest_ind(dm, nm)[1])

C_NM, C_DM = "#4C72B0", "#C0504D"
fig, axes = S.grid(1, 2, w=6.3, h=2.72, wspace=0.42,
                   left=0.115, right=0.975, top=0.855, bottom=0.155,
                   widths=[0.92, 1.20])

# -------------------------------------------------------------- (a) 散点
ax = axes[0]
rng = np.random.default_rng(7)
for i, (vals, col, lab) in enumerate([(nm, C_NM, "正常"), (dm, C_DM, "DM")]):
    jx = i + np.array([-0.085, 0.0, 0.085])
    ax.scatter(jx, vals, s=17, facecolor=col, edgecolor="white", lw=0.5,
               zorder=4, alpha=0.95)
    m, sd = vals.mean(), vals.std(ddof=1)
    ax.plot([i - 0.22, i + 0.22], [m, m], color="#222222", lw=1.25, zorder=5)
    ax.errorbar([i], [m], yerr=[sd], fmt="none", ecolor="#222222",
                elinewidth=0.65, capsize=2.0, capthick=0.65, zorder=3)

ax.plot([0, 0, 1, 1], [0.655, 0.700, 0.700, 0.655], color="#333333", lw=0.6)
ax.text(0.5, 0.712, "P = %.3f" % p_w, ha="center", va="bottom", fontsize=6.0)
ax.text(0.5, 0.995, "Welch t 检验, n = 3 / 3", ha="center", va="top",
        transform=ax.transAxes, fontsize=5.6, color="#555555")
ax.set_xlim(-0.45, 1.45)
ax.set_ylim(-0.62, 0.92)
ax.set_xticks([0, 1])
ax.set_xticklabels(["正常", "DM"], fontsize=6.8)
ax.set_yticks([-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6])
ax.set_ylabel("程序得分")
ax.set_title("上皮分化程序得分", pad=6)
ax.text(0.5, 0.03, "DM 组均值高 0.515", transform=ax.transAxes, ha="center",
        va="bottom", fontsize=5.8, color="#555555")
S.ygrid(ax)
S.clean(ax)
S.panel(ax, "a", x=-0.28, y=1.02)

# -------------------------------------------------------------- (b) dumbbell
ax = axes[1]
n_sorted = np.sort(nm)
d_sorted = np.sort(dm)
ypos = np.arange(len(d_sorted))[::-1]
for y, a, b in zip(ypos, n_sorted, d_sorted):
    up = b > a
    ax.annotate("", xy=(1.0, y), xytext=(0.0, y),
                arrowprops=dict(arrowstyle="-|>", lw=0.9, color=C_DM,
                                shrinkA=2.2, shrinkB=2.2, mutation_scale=7))
    ax.scatter([0], [y], s=17, facecolor=C_NM, edgecolor="white", lw=0.5, zorder=4)
    ax.scatter([1], [y], s=17, facecolor=C_DM, edgecolor="white", lw=0.5, zorder=4)
    ax.text(1.07, y, "%+.2f" % (b - a), va="center", ha="left", fontsize=6.0,
            color=C_DM if up else C_NM)
    ax.text(-0.07, y, "%.2f" % a, va="center", ha="right", fontsize=5.8,
            color="#555555")
ax.set_xlim(-0.30, 1.30)
ax.set_ylim(-0.55, len(ypos) - 0.45)
ax.set_xticks([0, 1])
ax.set_xticklabels(["正常", "DM"], fontsize=6.8)
ax.set_yticks(ypos)
ax.set_yticklabels(["排名 %d" % (i + 1) for i in range(len(ypos))], fontsize=6.0)
ax.set_ylabel("配对样本（按得分排名）")
ax.set_title("逐样本配对变化", pad=6)
S.clean(ax)
S.panel(ax, "b", x=-0.30, y=1.02)

S.save(fig, "S7_GSE326016")
print("deltas", np.round(d_sorted - n_sorted, 3))