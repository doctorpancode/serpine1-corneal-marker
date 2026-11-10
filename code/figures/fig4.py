# -*- coding: utf-8 -*-
"""图4 空间定位：AUCell 活性、基因-细胞类型相关与解卷积权重"""
import numpy as np
import pandas as pd
import figstyle as S

CN = {"Endothelium": "内皮", "Endothelium_Lymphatic": "淋巴管内皮",
      "Epithelium_Basal": "上皮基底", "Epithelium_Proliferating": "上皮增殖",
      "Epithelium_Superficial": "上皮浅层", "Epithelium_Wing": "上皮翼",
      "Fibroblast": "成纤维", "Macrophage": "巨噬",
      "Neutrophil": "中性粒", "Pericyte_SMC": "周细胞/SMC",
      "Schwann": "施万细胞"}
GROUPS = ["Control", "Keratoconus", "Viral", "Bacterial_fungal"]
GLAB = ["正常", "圆锥", "病毒", "细菌/真菌"]

au = S.load(6, "D6_aucell_summary.tsv")
cor = S.load(7, "D7_coregene_vs_celltype_cor.tsv")
dec = S.load(7, "D7_deconvolution_by_group.tsv")

fig, axes = S.grid(2, 2, w=6.5, h=4.30, wspace=0.62, hspace=0.62,
                   left=0.135, right=0.975, top=0.935, bottom=0.085)


def rho_panel(ax, gene, title, note=None):
    d = cor[cor.gene == gene].sort_values("rho")
    y = np.arange(len(d))
    cols = ["#C0504D" if v > 0 else "#5B8DB8" for v in d.rho]
    ax.axvline(0, color="#999999", lw=0.6)
    ax.hlines(y, 0, d.rho.values, color=cols, lw=1.1)
    ax.scatter(d.rho.values, y, s=13, c=cols, lw=0, zorder=3)
    for yi, (v, f) in enumerate(zip(d.rho, d.FDR)):
        ax.text(v + (0.05 if v > 0 else -0.05), yi, f"{v:+.2f}", va="center",
                ha="left" if v > 0 else "right", fontsize=5.2)
        if f < 0.05:
            ax.scatter([0.0], [yi], s=0, marker="")
    ax.set_yticks(y)
    ax.set_yticklabels([CN[c] for c in d.celltype], fontsize=5.8)
    ax.set_xlim(-1.0, 1.0)
    ax.set_xticks([-1, -0.5, 0, 0.5, 1])
    ax.set_xlabel("Spearman rho（21 例）")
    ax.set_title(title, pad=5)
    S.xgrid(ax); S.clean(ax)
    return d


# ------------------------------------------------- (a) SERPINE1 AUCell
ax = axes[0]
d = au[au.gene == "Serpine1"].pivot_table(index="celltype", columns="sample",
                                          values="mean_auc")
d = d.sort_values("NaOH")
y = np.arange(len(d))
h = 0.36
ax.barh(y + h / 2, d["uninjury"], height=h, color="#5B8DB8",
        edgecolor="none", label="未损伤")
ax.barh(y - h / 2, d["NaOH"], height=h, color="#C0504D",
        edgecolor="none", label="碱烧伤")
for yi, (a, b) in enumerate(zip(d["uninjury"], d["NaOH"])):
    ax.text(max(a, b) + 0.006, yi, f"{b:.3f} / {a:.3f}", va="center",
            fontsize=5.0, color="#333333")
ax.set_yticks(y); ax.set_yticklabels([CN[c] for c in d.index], fontsize=5.8)
ax.set_xlim(0, 0.345)
ax.set_xlabel("SERPINE1 平均 AUCell 活性")
ax.legend(loc="lower right", fontsize=5.6, frameon=False, handlelength=1.0,
          borderpad=0.2, labelspacing=0.25)
ax.annotate("59.7×", xy=(0.2234, len(d) - 1.30), xytext=(0.268, len(d) - 1.95),
            fontsize=5.8, color="#8C2F2C",
            arrowprops=dict(arrowstyle="-|>", lw=0.6, color="#8C2F2C"))
ax.set_title("SERPINE1 的细胞类型定位", pad=5)
S.xgrid(ax); S.clean(ax); S.panel(ax, "a", x=-0.42, y=1.02)

# ------------------------------------------------- (b) SERPINE1 rho
rho_panel(axes[1], "SERPINE1", "SERPINE1 与细胞丰度的关联")
S.panel(axes[1], "b", x=-0.34, y=1.02)

# ------------------------------------------------- (c) TMSB4X rho
rho_panel(axes[2], "TMSB4X", "TMSB4X 与细胞丰度的关联")
S.panel(axes[2], "c", x=-0.42, y=1.02)

# ------------------------------------------------- (d) deconvolution
ax = axes[3]
p = dec.pivot_table(index="celltype", columns="group", values="mean_w")
p = p.loc[sorted(p.index, key=lambda c: -p.loc[c].mean())][GROUPS]
S.heat(ax, p.values * 100, [CN[c] for c in p.index], GLAB,
       annot=np.round(p.values * 100, 1), annot_size=5.0,
       cmap=S.CMAP_SEQ, center=None, vmin=0, vmax=18.5,
       cbar_label="NNLS 权重 (%)", xrot=0,
       cbar_kw=dict(fraction=0.05, pad=0.03))
ax.set_xticks(range(4))
ax.set_xticklabels(GLAB, rotation=0, ha="center", fontsize=5.8)
ax.set_title("解卷积估计的细胞组成", pad=5)
S.panel(ax, "d", x=-0.38, y=1.02)

S.save(fig, "Fig4_localization")