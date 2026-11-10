# -*- coding: utf-8 -*-
"""图2 bulk 四组表达梯度、核心基因模式与差异表达
注：matplotlib 含 $...$ 的字符串整体交给 mathtext，中文字形会缺字，
   故含数学符号的文本一律只用 ASCII，中文文本一律不含 $。
"""
import numpy as np
import figstyle as S

fig, axes = S.grid(1, 3, w=6.4, h=2.72, wspace=0.62, left=0.085, right=0.975,
                   top=0.85, bottom=0.185, widths=[1.0, 1.12, 1.46])

groups = ["Control", "Keratoconus", "Viral", "Bacterial_fungal"]
lab = ["正常", "圆锥", "病毒性", "细菌/\n真菌"]
gs = S.load(5, "D5_group_stats.tsv")

# -------------------------------------------------- (a) SERPINE1 gradient
ax = axes[0]
sp = gs[gs.gene == "SERPINE1"].set_index("group").loc[groups]
x = np.arange(4)
ax.bar(x, sp["mean"], yerr=sp["sd"], color=S.BULK_COLORS, width=0.66,
       edgecolor="none", error_kw=dict(lw=0.6, capsize=2.0, ecolor="#555555"))
for xi, (m, s) in enumerate(zip(sp["mean"], sp["sd"])):
    ax.text(xi, m + s + 0.33, f"{m:.2f}", ha="center", fontsize=6.2)
ax.set_xticks(x)
ax.set_xticklabels([f"{l}\nn={int(v)}" for l, v in zip(lab, sp["n"])], fontsize=6.0)
ax.tick_params(axis="x", pad=2)
ax.set_ylabel("SERPINE1 表达量")
ax.set_ylim(0, 13.4)
ax.set_yticks([0, 3, 6, 9, 12])
ax.annotate("", xy=(3, 12.05), xytext=(0, 12.05),
            arrowprops=dict(arrowstyle="-|>", lw=0.7, color="#8A97A4",
                            shrinkA=0, shrinkB=0))
ax.text(1.5, 12.30, "5.4×   $P=9.4\\times10^{-4}$", ha="center", fontsize=6.0)
ax.set_title("SERPINE1 随病损程度递增", pad=6)
S.ygrid(ax); S.clean(ax); S.panel(ax, "a", x=-0.34)

# -------------------------------------------------- (b) core gene heatmap
ax = axes[1]
core = ["SERPINE1", "P4HB", "S100A13", "XIAP", "CAPG", "PRDX2",
        "VDAC1", "SOD1", "GJA1", "TMSB4X"]
mat = np.zeros((len(core), 4))
raw = np.zeros_like(mat)
for i, g in enumerate(core):
    v = gs[gs.gene == g].set_index("group").loc[groups]["mean"].values
    raw[i] = v
    mat[i] = (v - v.mean()) / v.std(ddof=0)
annot = [[f"{raw[i, j]:.1f}" for j in range(4)] for i in range(len(core))]
cols = ["正常", "圆锥", "病毒", "细菌"]
S.heat(ax, mat, core, cols, annot=annot, annot_size=5.6,
       cbar_label="行 z-score", xrot=0, cbar_kw=dict(fraction=0.052, pad=0.03))
ax.set_xticks(range(4))
ax.set_xticklabels(cols, rotation=0, ha="center", fontsize=6.0)
ax.set_title("10 个核心基因的组间模式", pad=6)
S.panel(ax, "b", x=-0.46)

# -------------------------------------------------- (c) volcano
ax = axes[2]
dg = S.load(3, "D3_limma_deg.tsv")
dg["nlp"] = -np.log10(dg["P"].clip(lower=1e-300))
sig = (dg.logFC.abs() > 0.5) & (dg.adjP < 0.05)
up = sig & (dg.logFC > 0)
dn = sig & (dg.logFC < 0)
n_up, n_dn = int(up.sum()), int(dn.sum())
print("DEG |logFC|>0.5 & adjP<0.05:", n_up + n_dn, "up", n_up, "down", n_dn)
ax.scatter(dg.loc[~sig, "logFC"], dg.loc[~sig, "nlp"], s=1.0, c="#DCDCDC",
           lw=0, rasterized=True)
ax.scatter(dg.loc[up, "logFC"], dg.loc[up, "nlp"], s=1.6, c=S.PAL["red"],
           lw=0, rasterized=True, label=f"上调 {n_up}")
ax.scatter(dg.loc[dn, "logFC"], dg.loc[dn, "nlp"], s=1.6, c=S.PAL["blue"],
           lw=0, rasterized=True, label=f"下调 {n_dn}")
# (gene, xoffset, yoffset, ha)
mark = [("SERPINE1", 7, -1, "left"), ("P4HB", 7, -9, "left"),
        ("XIAP", 8, 4, "left"), ("CAPG", -8, -7, "right"),
        ("S100A13", 8, 0, "left"), ("TMSB4X", 0, 9, "center"),
        ("GJA1", 7, -2, "left")]
for g, dx, dy, ha in mark:
    r = dg[dg.gene == g]
    if r.empty:
        continue
    r = r.iloc[0]
    ax.scatter([r.logFC], [r.nlp], s=10, facecolor="none",
               edgecolor="#222222", lw=0.6, zorder=5)
    ax.annotate(g, (r.logFC, r.nlp), textcoords="offset points",
                xytext=(dx, dy), fontsize=5.6, color="#111111", ha=ha)
for xv in (-0.5, 0.5):
    ax.axvline(xv, color="#999999", lw=0.45, ls=(0, (3, 2)))
ax.axhline(-np.log10(0.05), color="#999999", lw=0.45, ls=(0, (3, 2)))
ax.set_xlabel("$log_2$ fold change")
ax.set_ylabel("$-log_{10} P$".replace("$", "$", 2))
ax.set_xlim(-4.4, 5.2)
ax.set_ylim(0, 9.6)
ax.legend(loc="upper left", fontsize=5.8, markerscale=3.4, handletextpad=0.3,
          borderpad=0.25, labelspacing=0.25)
ax.text(0.985, 0.02, f"共 {n_up + n_dn} 个差异基因", transform=ax.transAxes,
        ha="right", va="bottom", fontsize=5.8, color="#555555")
ax.set_title("损伤 vs 正常的差异表达", pad=6)
S.clean(ax); S.panel(ax, "c", x=-0.24)

S.save(fig, "Fig2_bulk")