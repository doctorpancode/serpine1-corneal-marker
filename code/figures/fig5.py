# -*- coding: utf-8 -*-
"""图5 免疫浸润、细胞通讯与通路富集"""
import numpy as np
import pandas as pd
import figstyle as S

im = S.load(5, "D5_immune_group_stats.tsv")
cp = S.load(9, "D9_cellchat_compare_pathways.tsv")
lr = S.load(9, "D9_cellchat_lr_to_vessel_real.tsv")
kegg = S.load(5, "D5_enrich_KEGG.tsv")
reac = S.load(5, "D5_enrich_Reactome.tsv")

fig, axes = S.grid(2, 2, w=6.6, h=4.55, wspace=0.78, hspace=0.62,
                   left=0.155, right=0.975, top=0.935, bottom=0.085)

# ------------------------------------------- (a) immune infiltration delta
ax = axes[0]
sig = im[im.adjP < 0.05].sort_values("delta")
y = np.arange(len(sig))
cols = ["#C0504D" if v > 0 else "#5B8DB8" for v in sig.delta]
ax.barh(y, sig.delta, color=cols, height=0.62, edgecolor="none")
for yi, (dv, p) in enumerate(zip(sig.delta, sig.adjP)):
    ax.text(dv + 0.006, yi, f"{dv:+.3f}", va="center", fontsize=5.2)
    ax.text(-0.012, yi, S.stars(p), va="center", ha="right", fontsize=5.2,
            color="#555555")
ax.set_yticks(y)
ax.set_yticklabels([c.replace("_", " ") for c in sig.celltype], fontsize=5.6)
ax.set_xlim(-0.02, 0.30)
ax.set_xlabel("浸润分数差（损伤 − 正常）")
ax.set_title("免疫浸润改变（adjP < 0.05，n = 10）", pad=5)
S.xgrid(ax); S.clean(ax); S.panel(ax, "a", x=-0.62, y=1.02)

# ------------------------------------------- (b) CellChat pathways
ax = axes[1]
pick = ["SPP1", "TENASCIN", "CCL", "CXCL", "PTN", "PLAU", "TGFb", "VEGF",
        "LAMININ", "COLLAGEN", "FN1", "EGF"]
d = cp[cp.pathway.isin(pick)].copy()
d["ord"] = d.pathway.map({p: i for i, p in enumerate(pick)})
d = d.sort_values("log2FC")
y = np.arange(len(d))
cols = ["#C0504D" if v > 0 else "#5B8DB8" for v in d.log2FC]
ax.barh(y, d.log2FC, color=cols, height=0.62, edgecolor="none")
for yi, (v, r) in enumerate(zip(d.log2FC, d.ratio)):
    ax.text(max(v, 0) + 0.25, yi, f"×{r:,.1f}".replace(",", " "),
            va="center", fontsize=5.0, color="#333333")
ax.set_yticks(y); ax.set_yticklabels(d.pathway, fontsize=5.8)
ax.axvline(0, color="#777777", lw=0.6)
ax.set_xlim(-1.0, 18.2)
ax.set_xlabel("通讯强度倍数（损伤/正常，log2）")
ax.set_title("细胞通讯通路改变", pad=5)
S.xgrid(ax); S.clean(ax); S.panel(ax, "b", x=-0.34, y=1.02)

# ------------------------------------------- (c) ligand-receptor pairs
ax = axes[2]
AB = {"Neutrophil": "Neu", "Macrophage": "Mac", "Endothelium": "Endo",
      "Endothelium_Lymphatic": "LEC", "Fibroblast": "Fib",
      "Epithelium_Basal": "Basal", "Pericyte_SMC": "Peri"}
d = lr.nlargest(12, "prob").iloc[::-1]
y = np.arange(len(d))
src_col = {"Neutrophil": "#B3487E", "Macrophage": "#3B3B6D",
           "Fibroblast": "#4E9E8F", "Endothelium": "#2F6FA8"}
cols = [src_col.get(s, "#7A7A7A") for s in d.source]
ax.barh(y, d.prob * 1e3, color=cols, height=0.64, edgecolor="none")
labs = [f"{AB.get(s, s)}→{AB.get(t, t)}  {l}" for s, t, l in
        zip(d.source, d.target, d.lr)]
ax.set_yticks(y); ax.set_yticklabels(labs, fontsize=5.2)
for yi, v in enumerate(d.prob * 1e3):
    ax.text(v + 0.012, yi, f"{v:.2f}", va="center", fontsize=5.0)
ax.set_xlim(0, 0.62)
ax.set_xlabel("通讯概率（×1000）")
ax.set_title("指向血管内皮的配体-受体对", pad=5)
S.xgrid(ax); S.clean(ax); S.panel(ax, "c", x=-0.80, y=1.02)

# ------------------------------------------- (d) enrichment
ax = axes[3]
kpick = ["VEGF signaling pathway", "PI3K-Akt signaling pathway",
         "Focal adhesion", "Proteoglycans in cancer",
         "Lipid and atherosclerosis", "EGFR tyrosine kinase inhibitor resistance"]
rpick = ["REACTOME_SIGNALING_BY_VEGF",
         "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
         "REACTOME_DETOXIFICATION_OF_REACTIVE_OXYGEN_SPECIES",
         "REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION",
         "REACTOME_INTRINSIC_PATHWAY_FOR_APOPTOSIS",
         "REACTOME_APOPTOSIS"]


TITLE = {"REACTOME_SIGNALING_BY_VEGF": "Signaling by VEGF",
         "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION": "ECM organization",
         "REACTOME_DETOXIFICATION_OF_REACTIVE_OXYGEN_SPECIES": "ROS detoxification",
         "REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION": "Platelet activation",
         "REACTOME_INTRINSIC_PATHWAY_FOR_APOPTOSIS": "Intrinsic apoptosis",
         "REACTOME_APOPTOSIS": "Apoptosis"}


def short(s):
    return TITLE.get(s, s.replace("REACTOME_", "").replace("_", " ").title())


rows = []
for _, r in kegg[kegg.Description.isin(kpick)].iterrows():
    rows.append(("KEGG", r.Description, r.FoldEnrichment, r.Count, r["p.adjust"]))
for _, r in reac[reac.Description.isin(rpick)].iterrows():
    rows.append(("R", short(r.Description), r.FoldEnrichment, r.Count,
                 r["p.adjust"]))
d = pd.DataFrame(rows, columns=["db", "term", "fe", "n", "padj"])
d = d.sort_values(["db", "fe"])
y = np.arange(len(d))
c = -np.log10(d.padj.clip(lower=1e-9))
sc = ax.scatter(d.fe, y, s=d.n * 3.4, c=c, cmap="YlOrRd", lw=0.4,
                edgecolor="#555555", vmin=4, vmax=9)
for yi, (fe, n) in enumerate(zip(d.fe, d.n)):
    ax.text(fe + 0.55, yi, f"n={n}", va="center", fontsize=5.0, color="#444444")
ax.set_yticks(y); ax.set_yticklabels(d.term.tolist(), fontsize=5.4)
ax.set_xlim(0, 16.5)
ax.set_xlabel("Fold enrichment")
ax.set_title("损伤相关通路富集", pad=5)
cb = fig.colorbar(sc, ax=ax, fraction=0.045, pad=0.02)
cb.set_label("$-log_{10}$ adjP", fontsize=6.2)
cb.ax.tick_params(labelsize=5.6, width=0.5, length=2)
cb.outline.set_linewidth(0.4)
S.xgrid(ax); S.clean(ax); S.panel(ax, "d", x=-0.86, y=1.02)

S.save(fig, "Fig5_communication")
