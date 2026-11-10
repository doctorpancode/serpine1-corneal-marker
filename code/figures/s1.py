# -*- coding: utf-8 -*-
"""图S1  GSE58291 三个核心基因在对照 / FK / BK 三组的表达（mean ± SD）
均值取自 D9_core_genes_GSE58291.tsv（表为准）；SD 由归一化矩阵逐样本值计算
（Control/BK 用 BK_vs_Control 文件，FK 用 FK_vs_Control 文件，与表中口径一致）。
注：含 $...$ 的字符串只用 ASCII，中文串不含 $。
"""
import os
import numpy as np
import pandas as pd
import figstyle as S

GENES = ["SERPINE1", "P4HB", "GJA1"]
GROUPS = ["Control", "FK", "BK"]
GC = ["#4C72B0", "#E0A16B", "#C0504D"]
XL = ["对照\nn=12", "FK\nn=8", "BK\nn=7"]
YTICKS = {"SERPINE1": [8, 9, 10, 11, 12],
          "P4HB": [12.0, 12.5, 13.0, 13.5],
          "GJA1": [10, 11, 12, 13]}


def group_values(genes):
    """按样本取三基因的归一化表达值。"""
    sg = S.load(9, "D9_GSE58291_sample_group.tsv")
    gmap = dict(zip(sg["sample"], sg["group"]))
    pr = pd.read_csv(os.path.join(S.RAW, "GSE58291_expr.tsv"), sep="\t",
                     usecols=["PROBE_ID", "SYMBOL"], dtype=str)
    files = [("GSE58291_BK_vs_Control_normalized.tsv", {"Control", "BK"}),
             ("GSE58291_FK_vs_Control_normalized.tsv", {"FK"})]
    out = {}
    for fn, keep in files:
        mx = pd.read_csv(os.path.join(S.RAW, fn), sep="\t")
        cols = [c for c in mx.columns if c != "ID_REF" and "Detection" not in c]
        for g in genes:
            ids = pr.loc[pr.SYMBOL == g, "PROBE_ID"].tolist()
            r = mx.loc[mx.ID_REF.isin(ids), cols]
            if r.empty:
                continue
            v = pd.to_numeric(r.iloc[0], errors="coerce")
            d = out.setdefault(g, {k: [] for k in GROUPS})
            for c in cols:
                grp = gmap.get(c)
                if grp in keep:
                    d[grp].append(float(v[c]))
    return out


def pmtext(p):
    """ASCII mathtext 科学计数（含 $）。"""
    if p is None or p <= 0:
        return "$0$"
    e = int(np.floor(np.log10(p)))
    m = p / 10 ** e
    return f"${m:.1f}\\times10^{{{e}}}$"


tab = S.load(9, "D9_core_genes_GSE58291.tsv").set_index("gene").loc[GENES]
samp = group_values(GENES)
BB = dict(boxstyle="square,pad=0.16", facecolor="white", edgecolor="none")

fig, axes = S.grid(1, 3, w=6.4, h=2.55, wspace=0.40,
                   left=0.080, right=0.985, top=0.845, bottom=0.175)

for ax, g, letter in zip(axes, GENES, "abc"):
    m = np.array([tab.loc[g, f"mean_{k}"] for k in GROUPS], float)
    sd = np.array([np.std(samp[g][k], ddof=1) for k in GROUPS], float)
    span = float(max(m + sd) - min(m - sd))
    base = float(min(m - sd) - 0.10 * span)
    top = float(max(m + sd) + 0.46 * span)
    x = np.arange(3)
    ax.bar(x, m - base, bottom=base, width=0.62, color=GC,
           edgecolor="none", zorder=2)
    ax.errorbar(x, m, yerr=sd, fmt="none", ecolor="#4D4D4D", elinewidth=0.65,
                capsize=2.2, capthick=0.65, zorder=3)
    for xi, (mi, si) in enumerate(zip(m, sd)):
        ax.text(xi, mi + si + 0.020 * span, f"{mi:.2f}", ha="center",
                va="bottom", fontsize=6.0)
    ax.set_xlim(-0.55, 2.55)
    ax.set_ylim(base, top)
    ax.set_yticks(YTICKS[g])
    ax.set_yticklabels([("%.1f" % v) if g == "P4HB" else ("%d" % v)
                        for v in YTICKS[g]])
    ax.set_xticks(x)
    ax.set_xticklabels(XL, fontsize=5.8)
    ax.tick_params(axis="x", pad=2.0)
    S.ygrid(ax)
    S.clean(ax)
    ax.set_title(g, pad=5, fontsize=8.2)

    fc = float(tab.loc[g, "BK_logFC"])
    p = float(tab.loc[g, "BK_adjP"])
    rho = float(tab.loc[g, "trend_rho"])
    ax.text(0.5, 0.998, "BK log2FC = %+.2f" % fc,
            transform=ax.transAxes, ha="center", va="top",
            fontsize=5.5, color="#222222", bbox=BB)
    ax.text(0.5, 0.942, "adjP " + pmtext(p), transform=ax.transAxes,
            ha="center", va="top", fontsize=5.5, color="#222222", bbox=BB)
    ax.text(0.5, 0.884, "趋势 ρ = %+.2f" % rho, transform=ax.transAxes,
            ha="center", va="top", fontsize=5.5, color="#444444", bbox=BB)
    S.panel(ax, letter, x=-0.21, y=1.00)

axes[0].set_ylabel("表达量（log2 VST）")
S.save(fig, "S1_external_core")
print(tab[["mean_Control", "mean_FK", "mean_BK", "BK_logFC", "BK_adjP",
           "trend_rho"]].round(3).to_string())
for g in GENES:
    print(g, {k: round(float(np.std(samp[g][k], ddof=1)), 3) for k in GROUPS})