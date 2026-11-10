# -*- coding: utf-8 -*-
"""图S8  虚拟敲除打分的排名重排：原始权重 vs 行 L1 归一化
a 11 个核心基因的 rank_raw → rank_r1 dumbbell（横轴=排名，1 在左；
  每基因一行，行距均等，名次数字与基因名分列于两端，彻底避免排名
  1/50/92/136/139/175 这类相邻名次在纵轴刻度上压字）
b 三种相似度度量的 Spearman ρ（D9_regularization_compare.tsv）
"""
import numpy as np
import pandas as pd
import figstyle as S

CORE = ["TMSB4X", "PRDX2", "CAPG", "VDAC1", "GJA1", "P4HB", "VIM",
        "S100A13", "SOD1", "SERPINE1", "XIAP"]
C_UP, C_DN, C_EQ = "#3B6BA5", "#C0504D", "#9A9A9A"

d = S.load(9, "D9_regularized_ko_scores.tsv")
d["key"] = d.ko_gene.str.upper()
sub = d[d.key.isin(CORE)].copy()
N = len(d)

cmp = S.load(9, "D9_regularization_compare.tsv").set_index("metric")

fig, axes = S.grid(1, 2, w=6.4, h=4.05, wspace=0.60,
                   left=0.105, right=0.975, top=0.905, bottom=0.115,
                   widths=[1.52, 1.0])

# ---------------------------------------------- (a) 排名 dumbbell（横向）
ax = axes[0]
sub = sub.sort_values("rank_raw").reset_index(drop=True)
raw = sub.rank_raw.to_numpy(float)
new = sub.rank_r1.to_numpy(float)
names = sub.key.tolist()
ypos = np.arange(len(sub))[::-1].astype(float)

X_SEP, X_NUM = 985.0, 1012.0          # 名次数字列的分隔线与起点

for yv, a, b in zip(ypos, raw, new):
    dlt = b - a
    col = C_EQ if abs(dlt) < 1 else (C_UP if dlt < 0 else C_DN)
    if abs(dlt) < 1:
        ax.scatter([a], [yv], s=17, facecolor=col, edgecolor="white",
                   lw=0.5, zorder=5)
    else:
        ax.annotate("", xy=(b, yv), xytext=(a, yv),
                    arrowprops=dict(arrowstyle="-|>", lw=0.9, color=col,
                                    shrinkA=3.6, shrinkB=3.6, mutation_scale=7,
                                    alpha=0.95), zorder=3)
        ax.scatter([a], [yv], s=17, facecolor="white", edgecolor=col,
                   lw=0.8, zorder=5)
        ax.scatter([b], [yv], s=17, facecolor=col, edgecolor="white",
                   lw=0.5, zorder=5)
    ax.text(X_NUM, yv, "%d \u2192 %d" % (a, b), va="center", ha="left",
            fontsize=5.9, color=col)

ax.axvline(X_SEP, color="#E6E6E6", lw=0.6, zorder=0)
ax.set_xlim(-105, 1290)
ax.set_ylim(-1.55, len(sub) - 0.25)

# 图例（真实标记 + 彩色文字，不依赖字体里的圆圈符号）
KY = -1.02
ax.scatter([8], [KY], s=17, facecolor="white", edgecolor="#777777",
           lw=0.7, zorder=5)
ax.text(42, KY, "原始名次", fontsize=5.2, color="#777777", va="center",
        ha="left")
ax.scatter([255], [KY], s=17, facecolor="#777777", edgecolor="white",
           lw=0.5, zorder=5)
ax.text(289, KY, "行 L1 名次", fontsize=5.2, color="#777777", va="center",
        ha="left")
ax.text(600, KY, "前移", fontsize=5.2, color=C_UP, va="center", ha="left")
ax.text(690, KY, "·", fontsize=5.2, color="#999999", va="center", ha="left")
ax.text(715, KY, "后退", fontsize=5.2, color=C_DN, va="center", ha="left")
ax.text(805, KY, "·", fontsize=5.2, color="#999999", va="center", ha="left")
ax.text(830, KY, "不变", fontsize=5.2, color=C_EQ, va="center", ha="left")
ax.set_yticks(ypos)
ax.set_yticklabels(names, fontsize=6.3)
ax.set_xticks([0, 200, 400, 600, 800])
ax.set_xticklabels(["0", "200", "400", "600", "800"], fontsize=6.2)
ax.tick_params(axis="y", length=2.5, pad=2)
ax.set_xlabel("排名（1 = 最高，共 %d 个基因）" % N)
ax.xaxis.set_label_coords(0.40, -0.095)
ax.set_title("核心基因排名重排", pad=7)
ax.set_axisbelow(True)
ax.xaxis.grid(True, color="#EFEFEF", lw=0.5)
ax.yaxis.grid(False)
S.clean(ax)
S.panel(ax, "a", x=-0.15, y=1.055)

# ------------------------------------------------------- (b) ρ 条形图
ax = axes[1]
rows = [("rho_outdegree_vs_mean_raw", "出度 vs 原始权重"),
        ("rho_outdegree_vs_mean_rowL1norm", "出度 vs 行 L1 权重"),
        ("rho_outdegree_vs_cos_to_TMSB4X", "出度 vs TMSB4X 余弦")]
y = np.arange(len(rows))[::-1]
vals = [float(cmp.loc[k, "rho"]) for k, _ in rows]
ax.barh(y, vals, height=0.50, color="#4A6F96", edgecolor="none", zorder=2)
for yi, (k, _), v in zip(y, rows, vals):
    p = float(cmp.loc[k, "p"])
    ax.text(v + 0.028, yi + 0.11, "%.3f" % v, va="center", ha="left",
            fontsize=6.2)
    if p <= 0:
        txt = "$P < 0.001$"
    else:
        e = int(np.floor(np.log10(p)))
        m = p / 10 ** e
        txt = "$P = %.1f\\times10^{%d}$" % (m, e)
    ax.text(v + 0.028, yi - 0.15, txt, va="center", ha="left", fontsize=5.2,
            color="#666666")
ax.set_yticks(y)
ax.set_yticklabels([lab for _, lab in rows], fontsize=6.3)
ax.set_xlim(0, 1.20)
ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
ax.set_xlabel("Spearman ρ（与出度）")
ax.set_title("排名一致性", pad=7)
ax.set_axisbelow(True)
ax.xaxis.grid(True, color="#EDEDED", lw=0.5)
ax.yaxis.grid(False)
S.clean(ax)
S.panel(ax, "b", x=-0.46, y=1.055)

S.save(fig, "S8_rank_shift")
print(sub[["key", "rank_raw", "rank_r1"]].to_string(index=False))