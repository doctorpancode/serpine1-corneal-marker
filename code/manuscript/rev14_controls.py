# -*- coding: utf-8 -*-
"""Four control analyses for draft 1.4:
A housekeeping panel / global-shift test
B expression-decile matching for knockout effects
C self-effect removal from the epithelial 18-gene score
D secreted vs transcription-factor classes (GO-based, local GMT)
Writes _rev14_controls.json and prints a compact report.
"""
import json
import os
import re
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats

BASE = r"D:\A生信分析\虚拟敲除"
OUT = {}
ORDER = ["Control", "Keratoconus", "Viral", "Bacterial_fungal"]
EPI = ["Krt12", "Aqp5", "Gsta4", "Upk3bl", "Upk1b", "Krt15", "Krt17", "Krt19",
       "Krt13", "Dsg1a", "Perp", "Tmprss11d", "Lypd3", "Calml3", "Gsta2", "Sfn",
       "Krt5", "Krt14"]
CORE = ["Tmsb4x", "Serpine1", "P4hb", "S100a13", "Vim", "Prdx2", "Capg",
        "Vdac1", "Xiap", "Sod1", "Gja1"]

# ---------------------------------------------------------------- A
expr = pd.read_csv(os.path.join(BASE, r"D1执行结果\第一天执行结果\data\processed\GSE241715_expression.tsv"),
                   sep="\t", index_col=0)
sheet = pd.read_csv(os.path.join(BASE, r"D1执行结果\第一天执行结果\data\processed\GSE241715_sample_sheet.tsv"),
                    sep="\t").set_index("sample")
grp = sheet["group"]
lesion = grp != "Control"

HK = ["ACTB", "GAPDH", "RPLP0", "B2M", "PPIA", "EEF1A1", "UBC", "TUBB", "RPS18", "PSMB6"]
HK = [g for g in HK if g in expr.index]
hk_rows = {}
for g in HK:
    v = pd.Series({s: expr.loc[g, s] for s in grp.index})
    m = [v[grp == k].mean() for k in ORDER]
    lfc = m[-1] - m[0]
    hk_rows[g] = {"means": [round(x, 3) for x in m], "log2FC_BF_vs_Control": round(lfc, 3),
                  "sd": [round(v[grp == k].std(), 3) for k in ORDER]}
# all-gene distribution of BF-vs-Control log2FC
allfc = expr.loc[:, lesion].mean(1) - expr.loc[:, ~lesion].mean(1)
hk_pct = {g: round(float((allfc.abs() < abs(hk_rows[g]["log2FC_BF_vs_Control"])).mean() * 100), 1)
          for g in HK}
OUT["A_housekeeping"] = {
    "genes": hk_rows,
    "abs_log2FC_percentile_among_all_genes": hk_pct,
    "all_gene_log2FC_median": round(float(allfc.median()), 3),
    "all_gene_log2FC_iqr": [round(float(allfc.quantile(.25)), 3), round(float(allfc.quantile(.75)), 3)],
    "n_genes": int(allfc.size),
    "library": {c: [round(float(pd.to_numeric(sheet.loc[grp == k, c], errors="coerce").mean()), 1)
                    for k in ORDER]
                for c in ["library_raw_sum", "detected_genes_cpm_gt1", "median_log2cpm"] if c in sheet.columns},
}
# TMSB4X / SERPINE1 for reference
for g in ["TMSB4X", "SERPINE1"]:
    v = pd.Series({s: expr.loc[g, s] for s in grp.index})
    OUT["A_housekeeping"].setdefault("reference", {})[g] = {
        "means": [round(v[grp == k].mean(), 3) for k in ORDER],
        "log2FC_BF_vs_Control": round(v[grp == ORDER[-1]].mean() - v[grp == ORDER[0]].mean(), 3),
        "abs_log2FC_percentile": round(float((allfc.abs() < abs(v[grp == ORDER[-1]].mean() - v[grp == ORDER[0]].mean())).mean() * 100), 1),
    }

# ---------------------------------------------------------------- B
sig = pd.read_csv(os.path.join(BASE, r"第9天执行结果\results\D9_deconv_signature.tsv"),
                  sep="\t").set_index("gene").mean(axis=1)
sc = pd.read_csv(os.path.join(BASE, r"第9天执行结果\results\D9_regularized_ko_scores.tsv"),
                 sep="\t").set_index("ko_gene")
df = pd.DataFrame({"expr": sig, "eff": sc["mean_raw"], "deg": sc["out_degree"],
                   "resid_deg": sc["resid_raw"]})
df["dec"] = pd.qcut(df["expr"], 10, labels=False, duplicates="drop")
rho = lambda a, b: dict(zip(["rho", "p"], [round(float(x), 4) for x in stats.spearmanr(a, b)]))
dec_rows = {}
for g in CORE:
    d0 = int(df.loc[g, "dec"])
    sub = df[df.dec == d0]
    dec_rows[g] = {
        "expr": round(float(df.loc[g, "expr"]), 3),
        "expr_rank": int((df.expr > df.loc[g, "expr"]).sum()) + 1,
        "decile": d0, "decile_n": int(len(sub)),
        "decile_median_eff": round(float(sub.eff.median()), 5),
        "eff": round(float(df.loc[g, "eff"]), 6),
        "percentile_within_decile": round(float((sub.eff < df.loc[g, "eff"]).mean() * 100), 1),
    }
# two-stage detrend: degree (loess-like, via residual already available) then expression
x = df["expr"].values
y = df["resid_deg"].values
coef = np.polyfit(np.log1p(x), y, 2)
df["resid_deg_expr"] = y - np.polyval(coef, np.log1p(x))
df["rank_resid_deg_expr"] = df["resid_deg_expr"].rank(ascending=False).astype(int)
OUT["B_expression"] = {
    "spearman_eff_vs_expr": rho(df.eff, df.expr),
    "spearman_eff_vs_degree": rho(df.eff, df.deg),
    "spearman_resid_deg_vs_expr": rho(df.resid_deg, df.expr),
    "spearman_resid_both_vs_expr": rho(df.resid_deg_expr, df.expr),
    "decile_medians": {int(k): [int(v1), round(float(v2), 5), round(float(v3), 3)]
                       for k, v1, v2, v3 in zip(df.groupby("dec").size().index,
                                                df.groupby("dec").size().values,
                                                df.groupby("dec").eff.median().values,
                                                df.groupby("dec").expr.median().values)},
    "core": dec_rows,
    "rank_resid_degree_then_expr": {g: int(df.loc[g, "rank_resid_deg_expr"]) for g in CORE},
    "top10_after_both_detrend": [[g, round(float(df.loc[g, "resid_deg_expr"]), 5)]
                                 for g in df.resid_deg_expr.nlargest(10).index],
}

# ---------------------------------------------------------------- C
spec = pd.read_csv(os.path.join(BASE, r"第8天执行结果\results\D8_knockdown_spectrum.tsv"),
                   sep="\t", usecols=["ko_gene", "affected_gene", "distance", "direction"])
spec = spec[spec.affected_gene.isin(EPI)].copy()
spec["delta"] = spec.distance.abs() * spec.direction
wide = spec.pivot_table(index="ko_gene", columns="affected_gene", values="delta", aggfunc="sum")
wide = wide.reindex(columns=EPI)
inc = wide.abs().mean(axis=1)
exc = {}
for ko, row in wide.iterrows():
    r = row.drop(labels=[ko], errors="ignore")
    exc[ko] = float(r.abs().mean())
exc = pd.Series(exc)
rk_inc = inc.rank(ascending=False)
rk_exc = exc.rank(ascending=False)
OUT["C_self_effect"] = {
    "n_ko": int(len(inc)),
    "epi_genes_in_network": int(wide.index.isin(EPI).sum()),
    "epi_self_rank_before": {g: int(rk_inc[g]) for g in EPI if g in rk_inc.index},
    "epi_self_rank_after": {g: int(rk_exc[g]) for g in EPI if g in rk_exc.index},
    "median_epi_rank_before": float(np.median([rk_inc[g] for g in EPI if g in rk_inc.index])),
    "median_epi_rank_after": float(np.median([rk_exc[g] for g in EPI if g in rk_exc.index])),
    "core": {g: {"eff_incl": round(float(inc[g]), 6), "rank_incl": int(rk_inc[g]),
                 "eff_excl": round(float(exc[g]), 6), "rank_excl": int(rk_exc[g])}
             for g in CORE if g in inc.index},
    "top10_excl": [[g, round(float(exc[g]), 5)] for g in exc.nlargest(10).index],
}

# ---------------------------------------------------------------- D
# GO mapping is unavailable offline (all local .gmt files are 1-byte placeholders,
# msigdbr/msigdbdf not installed, org.Mm.eg.db absent). Use exhaustive symbol-family
# rules instead: applied to every gene in the network, not to a hand-picked subset.
SECRETED_RULES = [
    ("^Ccl", "chemokine"), ("^Cxcl", "chemokine"), ("^Cx3cl", "chemokine"),
    ("^Il[0-9]", "interleukin"), ("^Il1[0-9]", "interleukin"), ("^Il2[0-9]", "interleukin"),
    ("^Tgfb", "growth factor"), ("^Vegf", "growth factor"), ("^Igf", "growth factor"),
    ("^Fgf", "growth factor"), ("^Pdgfb?", "growth factor"), ("^Bmp", "growth factor"),
    ("^Gdf", "growth factor"), ("^Inhba", "growth factor"), ("^Lif$", "cytokine"),
    ("^Osm$", "cytokine"), ("^Mif$", "cytokine"), ("^Lcn2$", "secreted"),
    ("^Mmp", "protease"), ("^Timp", "protease inhibitor"), ("^Serpin", "protease inhibitor"),
    ("^Cst", "protease inhibitor"), ("^S100a(8|9|12)$", "alarmin"),
    ("^Apoe$", "secreted"), ("^Clu$", "secreted"), ("^Sparc", "matrix"),
    ("^Thbs", "matrix"), ("^Fn1$", "matrix"), ("^Col[0-9]", "matrix"),
    ("^Lama", "matrix"), ("^Lamb", "matrix"), ("^Lamc", "matrix"),
    ("^C3$", "complement"), ("^Cfb$", "complement"), ("^Cfh$", "complement"),
    ("^Pdgfa$", "growth factor"), ("^Wnt", "growth factor"),
]
TF_RULES = [
    ("^Fos", "AP-1"), ("^Jun", "AP-1"), ("^Atf", "AP-1"), ("^Egr", "immediate early"),
    ("^Klf", "Kruppel"), ("^Stat", "STAT"), ("^Nfkb", "NF-kB"), ("^Nfkbiz$", "NF-kB"),
    ("^Rela$", "NF-kB"), ("^Relb$", "NF-kB"), ("^Cebpb$", "C/EBP"), ("^Cebpd$", "C/EBP"),
    ("^Pax6$", "homeobox"), ("^Trp63$", "p53 family"), ("^Tp53$", "p53 family"),
    ("^Sox", "HMG"), ("^Snai", "zinc finger"), ("^Twist", "bHLH"), ("^Zeb", "zinc finger"),
    ("^Id[0-9]$", "bHLH"), ("^Hes", "bHLH"), ("^Myc$", "bHLH"), ("^Mycn$", "bHLH"),
    ("^Runx", "Runt"), ("^Smad", "SMAD"), ("^Nr4a", "nuclear receptor"),
    ("^Hif1a$", "bHLH-PAS"), ("^Epas1$", "bHLH-PAS"), ("^Yy1$", "zinc finger"),
    ("^Sp[0-9]$", "zinc finger"), ("^Ets", "ETS"), ("^Elk", "ETS"),
    ("^Xbp1$", "bZIP"), ("^Maf", "bZIP"), ("^Nfe2l2$", "bZIP"),
]


def classify(genes):
    out = {}
    for g in genes:
        s = next((f for p, f in SECRETED_RULES if re.match(p, g)), None)
        t = next((f for p, f in TF_RULES if re.match(p, g)), None)
        out[g] = ("secreted", s) if s and not t else ("tf", t) if t and not s else \
                 ("both", None) if s and t else ("neither", None)
    return out


net = set(sc.index)
cls = classify(net)
S = sorted(g for g in net if cls[g][0] == "secreted")
T = sorted(g for g in net if cls[g][0] == "tf")
N = sorted(g for g in net if cls[g][0] == "neither")
groups = {"secreted": S, "tf": T, "neither": N}
class_stats = {}
for k, gl in groups.items():
    v = sc.loc[gl, "mean_raw"]
    class_stats[k] = {"n": len(gl), "median_eff": round(float(v.median()), 6),
                      "median_rank_raw": float(sc.loc[gl, "rank_raw"].median()),
                      "median_expr": round(float(df.loc[gl, "expr"].median()), 3),
                      "median_eff_after_expr_detrend":
                          round(float(df.loc[gl, "resid_deg_expr"].median()), 6)}
# expression-matched null: compare each class against random gene sets of the same size
rng = np.random.default_rng(20261111)
ctrl = df["resid_deg_expr"]


def empirical_p(members, n_iter=5000):
    obs = float(df.loc[members, "resid_deg_expr"].median())
    vals = ctrl.values
    md = np.empty(n_iter)
    for i in range(n_iter):
        md[i] = np.median(rng.permutation(vals)[:len(members)])
    return round(obs, 6), round(float((md >= obs).mean()), 4), round(float(np.median(md)), 6)


OUT["D_classes"] = {
    "method": "symbol-family rules applied exhaustively to all 947 network genes (no GO mapping offline)",
    "class_stats": class_stats,
    "mannwhitney_secreted_lt_tf": dict(zip(["U", "p"], [float(x) for x in
                                                        stats.mannwhitneyu(sc.loc[S, "mean_raw"],
                                                                           sc.loc[T, "mean_raw"],
                                                                           alternative="less")])),
    "mannwhitney_after_expr_detrend": dict(zip(["U", "p"], [float(x) for x in
                                                            stats.mannwhitneyu(
                                                                df.loc[S, "resid_deg_expr"],
                                                                df.loc[T, "resid_deg_expr"],
                                                                alternative="less")])),
    "median_after_expr_detrend_secreted_vs_tf": [
        round(float(df.loc[S, "resid_deg_expr"].median()), 6),
        round(float(df.loc[T, "resid_deg_expr"].median()), 6)],
    "empirical_p_secreted": empirical_p(S),
    "empirical_p_tf": empirical_p(T),
    "secreted_members": [[g, round(float(sc.loc[g, "mean_raw"]), 6), int(sc.loc[g, "rank_raw"])] for g in S],
    "tf_members": [[g, round(float(sc.loc[g, "mean_raw"]), 6), int(sc.loc[g, "rank_raw"])] for g in T],
    "core_class": {g: ("secreted" if g in S else "tf" if g in T else "neither") for g in CORE},
}

with open(os.path.join(BASE, "_rev14_controls.json"), "w", encoding="utf-8") as fh:
    json.dump(OUT, fh, ensure_ascii=False, indent=1)

# ---------------------------------------------------------------- report
print("=== A housekeeping ===")
for g, d in OUT["A_housekeeping"]["genes"].items():
    print(f"{g:8s} means={d['means']} log2FC={d['log2FC_BF_vs_Control']:+.2f} pct={OUT['A_housekeeping']['abs_log2FC_percentile_among_all_genes'][g]}")
print("all-gene median log2FC", OUT["A_housekeeping"]["all_gene_log2FC_median"],
      "IQR", OUT["A_housekeeping"]["all_gene_log2FC_iqr"])
print("ref", json.dumps(OUT["A_housekeeping"]["reference"], ensure_ascii=False))
print("library", OUT["A_housekeeping"]["library"])

print("\n=== B expression deciles ===")
for k, v in OUT["B_expression"].items():
    if k in ("core", "top10_after_both_detrend", "decile_medians", "rank_resid_degree_then_expr"):
        continue
    print(k, v)
for g, d in OUT["B_expression"]["core"].items():
    print(f"{g:9s} expr={d['expr']:6.2f} rank={d['expr_rank']:4d} dec={d['decile']} decMed={d['decile_median_eff']:.5f} eff={d['eff']:.6f} pct={d['percentile_within_decile']}")
print("top10 after deg+expr detrend", OUT["B_expression"]["top10_after_both_detrend"])
print("core rank after both detrend", OUT["B_expression"]["rank_resid_degree_then_expr"])

print("\n=== C self-effect ===")
print("n_ko", OUT["C_self_effect"]["n_ko"], "epi in net", OUT["C_self_effect"]["epi_genes_in_network"])
print("median epi rank incl/excl", OUT["C_self_effect"]["median_epi_rank_before"], OUT["C_self_effect"]["median_epi_rank_after"])
print("top10 excl", OUT["C_self_effect"]["top10_excl"])
for g, d in OUT["C_self_effect"]["core"].items():
    print(f"{g:9s} incl={d['eff_incl']:.6f}/{d['rank_incl']:4d} excl={d['eff_excl']:.6f}/{d['rank_excl']:4d}")

print("\n=== C epi self ranks ===")
print("before", OUT["C_self_effect"]["epi_self_rank_before"])
print("after ", OUT["C_self_effect"]["epi_self_rank_after"])

print("\n=== D secreted vs TF ===")
print(json.dumps(OUT["D_classes"]["class_stats"], ensure_ascii=False))
print("MW", OUT["D_classes"]["mannwhitney_secreted_lt_tf"])
print("n secreted", len(S), "n tf", len(T), "n neither", len(N))
print("core class", OUT["D_classes"]["core_class"])
print("secreted sample", OUT["D_classes"]["secreted_members"][:12])
print("tf sample", OUT["D_classes"]["tf_members"][:12])