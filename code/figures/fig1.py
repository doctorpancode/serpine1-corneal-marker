# -*- coding: utf-8 -*-
"""图1 研究设计与分析流程"""
import numpy as np
import figstyle as S
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

W, H = 6.3, 5.05
fig = plt.figure(figsize=(W, H))
axt = fig.add_axes([0.01, 0.60, 0.98, 0.375])
axb = fig.add_axes([0.115, 0.075, 0.30, 0.40])
axc = fig.add_axes([0.60, 0.075, 0.365, 0.40])

# ------------------------------------------------ (a) pipeline
axt.set_xlim(0, 1); axt.set_ylim(0, 1); axt.axis("off")
S.panel(axt, "a", x=-0.005, y=0.94, size=11)

STAGE = ["数据输入", "网络与模块", "候选池与投票", "空间定位", "虚拟敲除", "拓扑校正与判读"]
COLX = [0.005, 0.175, 0.345, 0.515, 0.685, 0.855]
BW = 0.145

BOX = [
    # (col, ycenter, height, text, facecolor)
    (0, 0.68, 0.30, "GSE275678 单细胞\n19 493 × 15 873", "#EAF0F7"),
    (0, 0.32, 0.26, "GSE241715 bulk\n4 组 21 例", "#EAF0F7"),
    (0, 0.045, 0.20, "独立验证\nGSE58291 等 3 套", "#F4F4F4"),
    (1, 0.68, 0.30, "干性评分\n1.5007 → 1.3336\nP = 0.0056", "#EDF3EE"),
    (1, 0.30, 0.26, "WGCNA（power 10）\n6 模块 1 413 基因", "#EDF3EE"),
    (2, 0.68, 0.30, "三池交集 194\n并集 1 407", "#FBF1E6"),
    (2, 0.32, 0.26, "分池 T1–T3 = 167\nDEG 1 350 / 1 367", "#FBF1E6"),
    (2, 0.045, 0.18, "三算法计票 → 11 核心", "#FBF1E6"),
    (3, 0.68, 0.30, "AUCell 活性定位\nSERPINE1 / TMSB4X", "#F1EDF6"),
    (3, 0.32, 0.26, "NNLS 解卷积\n细胞权重 ρ", "#F1EDF6"),
    (3, 0.045, 0.18, "CellChat 2 029 互作", "#F1EDF6"),
    (4, 0.68, 0.30, "scTenifoldNet\n947 × 3 587\n边 600 604", "#FBEBEB"),
    (4, 0.30, 0.26, "scTenifoldKnk\n947 次敲除", "#FBEBEB"),
    (5, 0.68, 0.30, "出度三分位 / loess 残差\n行 L1 归一", "#FDF3E3"),
    (5, 0.36, 0.24, "上皮程序 mean|Δ|\n18 个基因", "#FDF3E3"),
    (5, 0.08, 0.18, "判读：标志 vs 节点", "#FDF3E3"),
]

for col, yc, hh, txt, fc in BOX:
    x = COLX[col]
    box = FancyBboxPatch((x, yc - hh / 2), BW, hh,
                         boxstyle="round,pad=0.004,rounding_size=0.012",
                         linewidth=0.6, edgecolor="#9AA5B1", facecolor=fc,
                         transform=axt.transAxes, clip_on=False)
    axt.add_patch(box)
    axt.text(x + BW / 2, yc, txt, ha="center", va="center", fontsize=5.0,
             linespacing=1.35, transform=axt.transAxes)

for i, s in enumerate(STAGE):
    axt.text(COLX[i] + BW / 2, 1.0, s, ha="center", va="bottom", fontsize=6.2,
             fontweight="bold", color="#31414F", transform=axt.transAxes)
    if i < len(STAGE) - 1:
        x0 = COLX[i] + BW + 0.004
        x1 = COLX[i + 1] - 0.004
        for yy in (0.62, 0.22):
            axt.add_patch(FancyArrowPatch((x0, yy), (x1, yy),
                                          arrowstyle="-|>", mutation_scale=6,
                                          lw=0.7, color="#8A97A4",
                                          transform=axt.transAxes, clip_on=False,
                                          shrinkA=0, shrinkB=0))

# ------------------------------------------------ (b) WGCNA module-trait heatmap
mt = S.load(2, "wgcna_module_trait_corr.tsv")
order = ["MEturquoise", "MEgreen", "MEred", "MEblue", "MEbrown", "MEyellow"]
mt = mt.set_index("module").loc[order]
mat = mt[["r_StemnessScore", "r_Lesion_binary"]].values
pv = mt[["p_StemnessScore", "p_Lesion_binary"]].values
annot = [[f"{mat[i, j]:.2f}{S.stars(pv[i, j])}" for j in range(2)] for i in range(len(order))]
S.heat(axb, mat, [f"{m.replace('ME','')} ({n})" for m, n in zip(order, mt.moduleSize)],
       ["干性评分", "病损"], annot=annot, annot_size=6.0,
       cbar_label="Pearson r", xrot=30)
axb.set_title("模块-性状关联", pad=6)
S.panel(axb, "b", x=-0.60, y=1.04)

# ------------------------------------------------ (c) three-pool venn counts
vn = S.load(3, "D3_venn_counts.tsv").set_index("region").n_genes
labels = ["CNV 相关集 A", "病理性血管\n新生集 B", "角膜疾病集 C",
          "A∩B", "A∩C", "B∩C", "A∩B∩C", "A∪B∪C"]
keys = ["A", "B", "C", "A&B", "A&C", "B&C", "A&B&C", "A|B|C"]
vals = [int(vn[k]) for k in keys]
cols = [S.PAL["blue"], S.PAL["teal"], S.PAL["amber"],
        "#9FB4C8", "#9FB4C8", "#9FB4C8", S.PAL["red"], "#B9B9B9"]
y = np.arange(len(labels))[::-1]
axc.barh(y, vals, color=cols, height=0.68, edgecolor="none")
axc.set_yticks(y); axc.set_yticklabels(labels, fontsize=6.4)
for yi, v in zip(y, vals):
    axc.text(v + 30, yi, f"{v}", va="center", fontsize=6.2)
axc.set_xlim(0, 1650)
axc.set_xlabel("基因数")
axc.set_title("候选基因池的交集结构", pad=6)
S.clean(axc); axc.xaxis.grid(True, color="#EDEDED", lw=0.5); axc.set_axisbelow(True)
S.panel(axc, "c", x=-0.34, y=1.02)

S.save(fig, "Fig1_design")


