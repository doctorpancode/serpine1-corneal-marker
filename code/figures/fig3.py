# -*- coding: utf-8 -*-
"""图3 单细胞图谱与损伤后的细胞组成重塑
a UMAP（11 类）；b 两组细胞占比配对条形；c NaOH/uninjured 占比对数倍数
"""
import numpy as np
import pandas as pd
import matplotlib.patheffects as pe
import figstyle as S

um = pd.read_csv(S.asset("D6_umap_meta.tsv"), sep="\t")
cw = S.load(6, "D6_celltype_composition_wide.tsv")

CN = {"Endothelium": "内皮", "Endothelium_Lymphatic": "淋巴管内皮",
      "Epithelium_Basal": "上皮基底", "Epithelium_Proliferating": "上皮增殖",
      "Epithelium_Superficial": "上皮浅层", "Epithelium_Wing": "上皮翼",
      "Fibroblast": "成纤维", "Macrophage": "巨噬",
      "Neutrophil": "中性粒", "Pericyte_SMC": "周细胞/SMC",
      "Schwann": "施万细胞"}
COL = {"Endothelium": "#2F6FA8", "Endothelium_Lymphatic": "#7B4EA3",
       "Epithelium_Basal": "#C0504D", "Epithelium_Proliferating": "#E08A3C",
       "Epithelium_Superficial": "#E8C547", "Epithelium_Wing": "#8FAF3F",
       "Fibroblast": "#4E9E8F", "Macrophage": "#3B3B6D",
       "Neutrophil": "#B3487E", "Pericyte_SMC": "#7A7A7A",
       "Schwann": "#5C3A21"}
order = ["Epithelium_Basal", "Epithelium_Wing", "Epithelium_Proliferating",
         "Epithelium_Superficial", "Fibroblast", "Neutrophil", "Macrophage",
         "Endothelium_Lymphatic", "Endothelium", "Pericyte_SMC", "Schwann"]

fig, axes = S.grid(1, 3, w=6.5, h=2.60, wspace=0.40, left=0.058, right=0.985,
                   top=0.86, bottom=0.15, widths=[1.28, 1.05, 0.85])

# ------------------------------------------------------------ (a) UMAP
ax = axes[0]
for ct in order:
    d = um[um.celltype == ct]
    ax.scatter(d.UMAP_1, d.UMAP_2, s=0.55, c=COL[ct], lw=0,
               rasterized=True, label=CN[ct])
OFF = {"Pericyte_SMC": (0.55, -0.10), "Endothelium": (-0.25, 0.10),
       "Macrophage": (-0.55, 0.35), "Neutrophil": (0.0, 0.25),
       "Endothelium_Lymphatic": (0.35, 0.30), "Schwann": (0.55, -0.20)}
for ct in order:
    d = um[um.celltype == ct]
    dx, dy = OFF.get(ct, (0.0, 0.0))
    cx, cy = d.UMAP_1.median() + dx, d.UMAP_2.median() + dy
    ax.text(cx, cy, CN[ct], fontsize=5.2, ha="center", va="center",
            color="#101010",
            path_effects=[pe.withStroke(linewidth=1.7, foreground="white")])
ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
ax.set_xticks([]); ax.set_yticks([])
ax.set_aspect("equal")
for s in ax.spines.values():
    s.set_visible(False)
ax.set_title("单细胞图谱（15 873 细胞）", pad=6)
S.panel(ax, "a", x=-0.04, y=1.10)

# ------------------------------------------------- (b) composition bars
ax = axes[1]
cw = cw.set_index("celltype").loc[order]
y = np.arange(len(order))
h = 0.36
ax.barh(y + h / 2, cw.frac_uninjury * 100, height=h, color="#5B8DB8",
        edgecolor="none", label="未损伤")
ax.barh(y - h / 2, cw.frac_NaOH * 100, height=h, color="#C0504D",
        edgecolor="none", label="碱烧伤")
for yi, (a, b) in enumerate(zip(cw.frac_uninjury * 100, cw.frac_NaOH * 100)):
    ax.text(max(a, b) + 0.9, yi, f"{b:.1f} / {a:.1f}", va="center",
            fontsize=5.3, color="#333333")
ax.set_yticks(y); ax.set_yticklabels([CN[c] for c in order], fontsize=6.0)
ax.invert_yaxis()
ax.set_xlim(0, 56)
ax.set_xlabel("细胞占比 (%)")
ax.legend(loc="lower right", fontsize=5.6, frameon=False, handlelength=1.1,
          borderpad=0.2, labelspacing=0.3)
ax.set_title("11 类细胞占比改变", pad=5)
S.xgrid(ax); S.clean(ax)
S.panel(ax, "b", x=-0.46)

# ------------------------------------------------- (c) log2 fold change
ax = axes[2]
eps = 1e-4
fc = np.log2((cw.frac_NaOH + eps) / (cw.frac_uninjury + eps))
ax.axvline(0, color="#999999", lw=0.6)
cols = ["#C0504D" if v > 0 else "#5B8DB8" for v in fc]
ax.hlines(y, 0, fc.values, color=cols, lw=1.1)
ax.scatter(fc.values, y, s=13, c=cols, zorder=3, lw=0)
for yi, v in enumerate(fc.values):
    t = "0" if abs(v) < 0.05 else f"{v:+.1f}"
    ax.text(v + (0.35 if v > 0 else -0.35), yi, t, va="center",
            ha="left" if v > 0 else "right", fontsize=5.3, color="#333333")
ax.set_yticks(y); ax.set_yticklabels(["" for _ in order])
ax.invert_yaxis()
ax.set_xlim(-9.9, 10.6)
ax.set_xlabel("倍数（碱烧伤/未损伤，log2）")
ax.set_title("组成变化幅度", pad=5)
S.xgrid(ax); S.clean(ax); ax.set_yticks([])
S.panel(ax, "c", x=-0.30)

S.save(fig, "Fig3_singlecell")
