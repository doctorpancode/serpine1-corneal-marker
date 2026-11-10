# -*- coding: utf-8 -*-
"""图S2  GSE58291 病损程度梯度：三基因在对照 / FK / BK 的均值斜率图
a 绝对值均值（mean ± SD，直接标注基因与趋势 ρ）
b 相对对照的变化量（等于表中的 FK_logFC / BK_logFC）
含 $...$ 串只用 ASCII，中文串不含 $。
"""
import os
import numpy as np
import pandas as pd
import figstyle as S

GENES = ["SERPINE1", "P4HB", "GJA1"]
GROUPS = ["Control", "FK", "BK"]
COL = {"SERPINE1": "#C0504D", "P4HB": "#E9A23B", "GJA1": "#3B6BA5"}
XL = ["对照", "FK", "BK"]
# 右侧直接标签的纵向微调（SERPINE1 与 GJA1 在 BK 端几乎重合）
DY = {"SERPINE1": -11, "P4HB": 0, "GJA1": 11}


def group_values(genes):
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


tab = S.load(9, "D9_core_genes_GSE58291.tsv").set_index("gene").loc[GENES]
samp = group_values(GENES)

fig, axes = S.grid(1, 2, w=6.4, h=2.85, wspace=0.42,
                   left=0.095, right=0.815, top=0.855, bottom=0.165)

# ------------------------------------------------------------ (a) 绝对值
ax = axes[0]
x = np.arange(3)
for g in GENES:
    m = np.array([tab.loc[g, "mean_%s" % k] for k in GROUPS], float)
    sd = np.array([np.std(samp[g][k], ddof=1) for k in GROUPS], float)
    ax.errorbar(x, m, yerr=sd, color=COL[g], lw=1.15, marker="o", ms=3.6,
                mfc="white", mew=1.0, ecolor=COL[g], elinewidth=0.6,
                capsize=1.8, zorder=3,
                label="%s   ρ = %+.2f" % (g, tab.loc[g, "trend_rho"]))
ax.axvline(0, color="#BBBBBB", lw=0.45, ls=(0, (2, 2)), zorder=1)
ax.axvline(1, color="#BBBBBB", lw=0.45, ls=(0, (2, 2)), zorder=1)
ax.set_xlim(-0.30, 2.30)
ax.set_ylim(7.28, 14.1)
ax.set_xticks(x)
ax.set_xticklabels(XL, fontsize=6.6)
ax.set_yticks([8, 9, 10, 11, 12, 13, 14])
ax.set_ylabel("表达量（log2 VST）")
ax.set_title("三组均值梯度（mean ± SD）", pad=6)
lg = ax.legend(loc="lower right", bbox_to_anchor=(1.005, -0.015), fontsize=5.8,
               labelspacing=0.28, handlelength=1.3, handletextpad=0.45,
               borderpad=0.1)
for t, g in zip(lg.get_texts(), GENES):
    t.set_color(COL[g])
S.ygrid(ax)
S.clean(ax)
S.panel(ax, "a", x=-0.13, y=1.02)

# ------------------------------------------------------------ (b) 变化量
ax = axes[1]
xg = np.array([0.0, 0.62, 1.24])
for i, g in enumerate(GENES):
    dv = np.array([0.0, tab.loc[g, "FK_logFC"], tab.loc[g, "BK_logFC"]], float)
    ax.plot(xg, dv, color=COL[g], lw=1.15, marker="o", ms=3.6, mfc="white",
            mew=1.0, zorder=3)
    ax.text(xg[2] + 0.06, dv[2], "%s\n%+.2f" % (g, dv[2]), color=COL[g],
            fontsize=5.7, va="center", ha="left", linespacing=1.35,
            clip_on=False)
ax.axhline(0, color="#8C8C8C", lw=0.6, zorder=1)
ax.set_xlim(-0.28, 1.24)
ax.set_ylim(-2.75, 3.55)
ax.set_xticks(xg)
ax.set_xticklabels(XL, fontsize=6.6)
ax.set_yticks([-2, -1, 0, 1, 2, 3])
ax.set_ylabel("相对对照变化（log2）")
ax.set_title("相对对照的变化方向", pad=6)
ax.text(0.5, 0.03, "SERPINE1 / P4HB 单调上升，GJA1 下降后再回升",
        transform=ax.transAxes, ha="center", va="bottom", fontsize=5.6,
        color="#555555")
S.ygrid(ax)
S.clean(ax)
S.panel(ax, "b", x=-0.16, y=1.02)

S.save(fig, "S2_gradient")
for g in GENES:
    print(g, "SD", {k: round(float(np.std(samp[g][k], ddof=1)), 3) for k in GROUPS},
          "FK_logFC", round(float(tab.loc[g, "FK_logFC"]), 3),
          "BK_logFC", round(float(tab.loc[g, "BK_logFC"]), 3))