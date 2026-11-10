# -*- coding: utf-8 -*-
"""图7 核效应量级：SERPINE1 作为标志而非节点"""
import numpy as np
import pandas as pd
import figstyle as S

mk = S.load(9, "D9_multiko_epi_table.tsv")
reg = S.load(9, "D9_regularized_ko_scores.tsv")
summ = S.load(9, "D9_multiko_summary.tsv")

KOS = ["SERPINE1", "P4HB", "TMSB4X"]
d = mk.copy()
d["mag"] = d[KOS].abs().max(axis=1)
d = d.sort_values("TMSB4X")

fig, axes = S.grid(2, 2, w=6.9, h=5.5, wspace=0.60, hspace=0.30,
                   left=0.165, right=0.975, top=0.925, bottom=0.095)

# ------------------------------------------------ (a) heatmap
ax = axes[0]
mat = d[KOS].values
ann = [[(f"{v:.1e}" if 0 < abs(v) < 1e-3 else f"{v:.3f}") for v in row]
       for row in mat]
S.heat(ax, mat, d.gene.tolist(), KOS, annot=ann,
       annot_size=5.6, cbar_label="上皮基因 Δ", xrot=0,
       cbar_kw=dict(fraction=0.05, pad=0.03))
ax.set_xticks(range(3))
ax.set_xticklabels(KOS, rotation=0, ha="center", fontsize=7.0)
ax.set_title("18 个上皮分化基因的敲除效应", pad=5)
S.panel(ax, "a", x=-0.50, y=1.02)

# ------------------------------------------------ (b) magnitude bars
ax = axes[1]
m = [summ.loc[summ.combo == k, "epi_mean_abs_multi"].iloc[0] for k in KOS]
x = np.arange(3)
ax.bar(x, m, color=["#C0504D", "#7E9DB8", "#3B3B6D"], width=0.55,
       edgecolor="none")
ax.set_yscale("log")
for xi, v in enumerate(m):
    ax.text(xi, v * 1.35, f"{v:.6f}", ha="center", fontsize=5.4)
ax.set_xticks(x); ax.set_xticklabels(KOS, fontsize=6.2)
ax.set_ylabel("上皮程序平均 |Δ|")
ax.set_ylim(2e-4, 1.6)
ax.annotate("", xy=(2, 0.30), xytext=(0, 0.30),
            arrowprops=dict(arrowstyle="<->", lw=0.6, color="#8A97A4"))
ax.text(1.0, 0.36, "480×", ha="center", fontsize=6.2, color="#333333")
ax.set_title("核效应的量级差", pad=5)
S.ygrid(ax); S.clean(ax); S.panel(ax, "b", x=-0.30, y=1.02)

# ------------------------------------------------ (c) per-gene scatter
ax = axes[2]
ax.scatter(d.SERPINE1, d.TMSB4X, s=11, c="#7E9DB8", lw=0.4,
           edgecolor="white", zorder=3)
lbl = ["Krt12", "Aqp5", "Upk3bl", "Krt5", "Perp", "Krt14"]
for g in lbl:
    r = d[d.gene == g]
    if r.empty:
        continue
    r = r.iloc[0]
    ax.annotate(g, (r.SERPINE1, r.TMSB4X), textcoords="offset points",
                xytext=(4, 3), fontsize=5.0)
ax.axhline(0, color="#BBBBBB", lw=0.5)
ax.axvline(0, color="#BBBBBB", lw=0.5)
ax.set_xlabel("SERPINE1 敲除的 Δ")
ax.set_ylabel("TMSB4X 敲除的 Δ")
ax.text(0.40, 0.06, "Pearson r = 0.096\n上皮内 Spearman = 0.569",
        transform=ax.transAxes, fontsize=5.4, va="bottom", color="#333333")
ax.set_title("两种敲除的效应方向不一致", pad=5)
S.clean(ax); S.panel(ax, "c", x=-0.30, y=1.02)

# ------------------------------------------------ (d) genome-wide context
ax = axes[3]
v = reg.mean_r1.clip(lower=1e-9)
ax.hist(np.log10(v), bins=48, color="#CFD6DE", edgecolor="white", lw=0.25)
for k, c in [("Serpine1", "#C0504D"), ("P4hb", "#7E9DB8"), ("Tmsb4x", "#3B3B6D")]:
    r = reg[reg.ko_gene == k]
    if r.empty:
        continue
    xv = np.log10(max(r.mean_r1.iloc[0], 1e-9))
    ax.axvline(xv, color=c, lw=0.9)
    ax.text(xv, ax.get_ylim()[1] * 0.92, k, rotation=90, fontsize=5.2,
            color=c, ha="right", va="top")
ax.set_xlabel("平均 |Δ|（log10，947 次敲除）")
ax.set_ylabel("敲除数")
ax.set_title("SERPINE1 位于效应分布的下端", pad=5)
S.ygrid(ax); S.clean(ax); S.panel(ax, "d", x=-0.26, y=1.02)

S.save(fig, "Fig7_nuclear_effect")
