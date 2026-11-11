# -*- coding: utf-8 -*-
"""Emit archive-ready tables for the four control analyses of draft 1.4."""
import io
import json
import os
import re

import numpy as np
import pandas as pd

BASE = r"D:\A生信分析\虚拟敲除"
D9 = os.path.join(BASE, r"第9天执行结果\results")

sc = pd.read_csv(os.path.join(D9, "D9_regularized_ko_scores.tsv"), sep="\t").set_index("ko_gene")
sig = pd.read_csv(os.path.join(D9, "D9_deconv_signature.tsv"), sep="\t").set_index("gene").mean(axis=1)

df = pd.DataFrame({"expr": sig}).join(sc, how="inner")
df["expr_rank"] = df["expr"].rank(ascending=False).astype(int)
df["decile"] = pd.qcut(df["expr"], 10, labels=False, duplicates="drop") + 1
df["pct_within_decile"] = df.groupby("decile")["mean_raw"].rank(pct=True) * 100
coef = np.polyfit(np.log1p(df["expr"].values), df["resid_raw"].values, 2)
df["resid_both"] = df["resid_raw"].values - np.polyval(coef, np.log1p(df["expr"].values))
df["rank_resid_both"] = df["resid_both"].rank(ascending=False).astype(int)

SECRETED_RULES = [
    ("^Ccl", "chemokine"), ("^Cxcl", "chemokine"), ("^Cx3cl", "chemokine"),
    ("^Il[0-9]", "interleukin"), ("^Tgfb", "growth factor"), ("^Vegf", "growth factor"),
    ("^Igf", "growth factor"), ("^Fgf", "growth factor"), ("^Pdgfb?", "growth factor"),
    ("^Bmp", "growth factor"), ("^Gdf", "growth factor"), ("^Inhba", "growth factor"),
    ("^Lif$", "cytokine"), ("^Osm$", "cytokine"), ("^Mif$", "cytokine"),
    ("^Lcn2$", "secreted"), ("^Mmp", "protease"), ("^Timp", "protease inhibitor"),
    ("^Serpin", "protease inhibitor"), ("^Cst", "protease inhibitor"),
    ("^S100a(8|9|12)$", "alarmin"), ("^Apoe$", "secreted"), ("^Clu$", "secreted"),
    ("^Sparc", "matrix"), ("^Thbs", "matrix"), ("^Fn1$", "matrix"), ("^Col[0-9]", "matrix"),
    ("^Lama", "matrix"), ("^Lamb", "matrix"), ("^Lamc", "matrix"),
    ("^C3$", "complement"), ("^Cfb$", "complement"), ("^Cfh$", "complement"),
    ("^Wnt", "growth factor"),
]
TF_RULES = [
    ("^Fos", "AP-1"), ("^Jun", "AP-1"), ("^Atf", "AP-1"), ("^Egr", "immediate early"),
    ("^Klf", "Kruppel"), ("^Stat", "STAT"), ("^Nfkb", "NF-kB"), ("^Rela$", "NF-kB"),
    ("^Relb$", "NF-kB"), ("^Cebpb$", "C/EBP"), ("^Cebpd$", "C/EBP"), ("^Pax6$", "homeobox"),
    ("^Trp63$", "p53 family"), ("^Tp53$", "p53 family"), ("^Sox", "HMG"),
    ("^Snai", "zinc finger"), ("^Twist", "bHLH"), ("^Zeb", "zinc finger"), ("^Id[0-9]$", "bHLH"),
    ("^Hes", "bHLH"), ("^Myc$", "bHLH"), ("^Mycn$", "bHLH"), ("^Runx", "Runt"),
    ("^Smad", "SMAD"), ("^Nr4a", "nuclear receptor"), ("^Hif1a$", "bHLH-PAS"),
    ("^Epas1$", "bHLH-PAS"), ("^Yy1$", "zinc finger"), ("^Sp[0-9]$", "zinc finger"),
    ("^Ets", "ETS"), ("^Elk", "ETS"), ("^Xbp1$", "bZIP"), ("^Maf", "bZIP"),
    ("^Nfe2l2$", "bZIP"),
]


def classify(g):
    s = next((f for p, f in SECRETED_RULES if re.match(p, g)), None)
    t = next((f for p, f in TF_RULES if re.match(p, g)), None)
    if s and not t:
        return "secreted", s
    if t and not s:
        return "tf", t
    if s and t:
        return "both", s
    return "neither", ""


cls = [classify(g) for g in df.index]
df["class"] = [c[0] for c in cls]
df["family"] = [c[1] for c in cls]

cols = ["expr", "expr_rank", "decile", "pct_within_decile", "out_degree", "mean_raw", "rank_raw",
        "mean_r1", "rank_r1", "cos_to_tmsb4x", "rank_cos", "resid_raw", "rank_resid",
        "resid_both", "rank_resid_both", "class", "family"]
out = df[cols].round(6).sort_values("rank_raw")
out.index.name = "ko_gene"
out.to_csv(os.path.join(D9, "D9_control_analyses.tsv"), sep="\t")
print("wrote D9_control_analyses.tsv", out.shape)
print(out["class"].value_counts().to_dict())
print("rules:", len(SECRETED_RULES), len(TF_RULES))

summary = {
    "n_genes": int(len(df)),
    "class_counts": df["class"].value_counts().to_dict(),
    "class_median_eff": df.groupby("class")["mean_raw"].median().round(6).to_dict(),
    "class_median_rank": df.groupby("class")["rank_raw"].median().to_dict(),
    "class_median_expr": df.groupby("class")["expr"].median().round(3).to_dict(),
    "spearman": {
        "eff_vs_expr": round(float(df["mean_raw"].corr(df["expr"], method="spearman")), 4),
        "eff_vs_degree": round(float(df["mean_raw"].corr(df["out_degree"], method="spearman")), 4),
        "resid_vs_expr": round(float(df["resid_raw"].corr(df["expr"], method="spearman")), 4),
        "resid_both_vs_expr": round(float(df["resid_both"].corr(df["expr"], method="spearman")), 4),
    },
    "core": {g: {k: (round(float(v), 6) if isinstance(v, (int, float, np.integer, np.floating))
                      else str(v))
                 for k, v in df.loc[g, cols].items()} for g in
             ["Tmsb4x", "Serpine1", "P4hb", "S100a13", "Vim", "Prdx2", "Capg", "Vdac1",
              "Xiap", "Sod1", "Gja1"]},
}
with io.open(os.path.join(BASE, "_rev14_tables_summary.json"), "w", encoding="utf-8") as fh:
    json.dump(summary, fh, ensure_ascii=False, indent=1)
print(json.dumps(summary["class_median_eff"], ensure_ascii=False),
      json.dumps(summary["spearman"], ensure_ascii=False))