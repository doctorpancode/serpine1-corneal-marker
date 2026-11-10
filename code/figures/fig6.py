# -*- coding: utf-8 -*-
"""图6 虚拟敲除的扰动幅度、网络位置与上皮选择性"""
import numpy as np
import pandas as pd
import figstyle as S

two = S.load(7, "D7_knk_official_twoKO.tsv")
top = S.load(8, "D8_top10_shared_targets.tsv")
sel = S.load(8, "D8_epi_selectivity_ratio.tsv")
reg = S.load(9, "D9_regularized_ko_scores.tsv")
res = S.load(9, "D9_core_genes_residual_rank.tsv")
strata = S.load(9, "D9_centrality_by_stratum.tsv")

CORE = ["Tmsb4x", "Prdx2", "Capg", "Vdac1", "Gja1", "P4hb", "Vim",
        "S100a13", "Sod1", "Serpine1", "Xiap"]

fig, axes = S.grid(2, 3, w=6.7, h=4.55, wspace=0.62, hspace=0.78,
                   left=0.088, right=0.985, top=0.935, bottom=0.095,
                   widths=[1.02, 1.0, 1.10])

# ------------------------------------------------ (a) two KO on epithelial
ax = axes[0]
genes = ["Krt12", "Aqp5", "Gsta4", "Upk3bl"]
d = two[two.gene.isin(genes)]
y = np.arange(len(genes))
h = 0.34
neg = d[d.target == "Serpine1"].set_index("gene").loc[genes]
tms = d[d.target == "Tmsb4x"].set_index("gene").loc[genes]
ax.barh(y + h / 2, neg.Z, height=h, color="#C0504D", edgecolor="none",
        label="Serpine1 KO")
ax.barh(y - h / 2, tms.Z, height=h, color="#2F6FA8", edgecolor="none",
        label="Tmsb4x KO")
for yi, g in enumerate(genes):
    ax.text(neg.Z[g] + 0.06, yi + h / 2, f"{neg.Z[g]:.2f}", va="center", fontsize=5.1)
    ax.text(tms.Z[g] + 0.06, yi - h / 2, f"{tms.Z[g]:.2f}", va="center", fontsize=5.1)
ax.set_yticks(y); ax.set_yticklabels(genes, fontsize=6.2)
ax.invert_yaxis()
ax.set_xlim(0, 4.6)
ax.set_xlabel("扰动强度 Z")
ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.005), fontsize=5.4,
          frameon=False, handlelength=1.0, borderpad=0.2, labelspacing=0.25,
          ncol=2)
ax.set_title("两个核心基因敲除的上皮响应", pad=16)
S.xgrid(ax); S.clean(ax); S.panel(ax, "a", x=-0.30, y=1.03)

# ------------------------------------------------ (b) shared targets
ax = axes[1]
d = top.nlargest(10, "N").sort_values("N")
y = np.arange(len(d))
ax.barh(y, d.N, color="#7E9DB8", height=0.62, edgecolor="none")
for yi, v in enumerate(d.N):
    ax.text(v + 5, yi, str(int(v)), va="center", fontsize=5.2)
ax.set_yticks(y); ax.set_yticklabels(d.affected_gene, fontsize=5.8)
ax.set_xlim(0, 500)
ax.set_xlabel("受扰动 KO 数")
ax.set_title("公共下游靶基因（前 10）", pad=5)
S.xgrid(ax); S.clean(ax); S.panel(ax, "b", x=-0.36, y=1.03)

# ------------------------------------------------ (c) degree vs effect
ax = axes[2]
d = reg.merge(sel[["ko_gene", "epithel"]], on="ko_gene", how="left")
d["epi"] = d["epithel"].clip(lower=1e-6)
ax.scatter(d.out_degree, d.epi, s=2.2, c="#C9CFD6", lw=0, rasterized=True)
hi = d[d.ko_gene.isin(CORE)]
ax.scatter(hi.out_degree, hi.epi, s=16, c="#C0504D", lw=0.3,
           edgecolor="white", zorder=4)
LB = {"Tmsb4x": (5, 2), "Prdx2": (5, 2), "Capg": (5, -6),
      "Serpine1": (5, -7), "Xiap": (5, 3)}
for _, r in hi.iterrows():
    if r.ko_gene in LB:
        ax.annotate(r.ko_gene, (r.out_degree, r.epi),
                    textcoords="offset points", xytext=LB[r.ko_gene],
                    fontsize=5.0, color="#111111")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("网络出度")
ax.set_ylabel("上皮 18 基因平均 |Δ|")
ax.text(0.04, 0.94, "Spearman rho = 0.917", transform=ax.transAxes,
        fontsize=5.8, va="top", color="#333333")
ax.set_title("扰动幅度取决于网络位置", pad=5)
S.clean(ax); S.panel(ax, "c", x=-0.32, y=1.03)

# ------------------------------------------------ (d) loess residual
ax = axes[3]
d = res.sort_values("epi_resid")
y = np.arange(len(d))
cols = ["#C0504D" if v > 0 else "#5B8DB8" for v in d.epi_resid]
ax.barh(y, d.epi_resid, color=cols, height=0.62, edgecolor="none")
for yi, (v, rk) in enumerate(zip(d.epi_resid, d.core_rank_resid)):
    lab = f"{v:+.4f}" if abs(v) < 0.01 else f"{v:+.3f}"
    ax.text(max(v, 0) + 0.004, yi, lab + f" (#{int(rk)})", va="center",
            ha="left", fontsize=4.9)
ax.set_yticks(y); ax.set_yticklabels(d.ko_gene, fontsize=5.8)
ax.axvline(0, color="#777777", lw=0.6)
ax.set_xlim(-0.008, 0.30)
ax.set_xlabel("残差（观测 − 期望）")
ax.set_title("剔除出度后的上皮效应残差", pad=5)
S.xgrid(ax); S.clean(ax); S.panel(ax, "d", x=-0.30, y=1.03)

# ------------------------------------------------ (e) selectivity test
ax = axes[4]
a = sel[sel.grp == "background"].ratio
b = sel[sel.grp == "core candidate"].ratio
bins = np.linspace(0, 3.0, 34)
ax.hist(a.clip(upper=3), bins=bins, color="#C9CFD6", edgecolor="white",
        lw=0.3, density=True, label="背景 KO（n = 936）")
ax.hist(b.clip(upper=3), bins=np.linspace(0, 3.0, 9), color="#C0504D",
        alpha=0.75, edgecolor="white", lw=0.3, density=True,
        label="候选基因（n = 11）")
ax.axvline(a.median(), color="#4A5A6A", lw=0.7, ls=(0, (3, 2)))
ax.axvline(b.median(), color="#8C2F2C", lw=0.9)
ax.text(0.03, 0.995, f"候选中位 {b.median():.2f} vs 背景 {a.median():.2f}\n"
        "Wilcoxon P = 0.676，置换 P = 0.413\nCliff δ = −0.080",
        transform=ax.transAxes, fontsize=5.3, va="top", color="#333333")
ax.set_xlabel("上皮/非上皮效应比")
ax.set_ylabel("相对频率")
ax.legend(loc="lower right", fontsize=5.2, frameon=False, handlelength=1.0,
          borderpad=0.15, labelspacing=0.25)
ax.set_title("上皮选择性检验", pad=5)
S.ygrid(ax); S.clean(ax); S.panel(ax, "e", x=-0.34, y=1.03)

# ------------------------------------------------ (f) stratified rho
ax = axes[5]
lab = {"low": "低出度层", "mid": "中出度层", "high": "高出度层"}
d = strata.copy()
d["ord"] = d.strata.map({"low": 0, "mid": 1, "high": 2})
d = d.sort_values("ord")
x = np.arange(len(d))
ax.bar(x, d.rho, color="#7E9DB8", width=0.55, edgecolor="none")
for xi, (v, n) in enumerate(zip(d.rho, d.n)):
    ax.text(xi, v + 0.018, f"{v:.3f}", ha="center", fontsize=5.6)
    ax.text(xi, 0.05, f"n={int(n)}", ha="center", fontsize=5.2, color="white")
ax.set_xticks(x); ax.set_xticklabels([lab[s] for s in d.strata], fontsize=6.0)
ax.set_ylim(0, 1.12)
ax.set_ylabel("Spearman rho")
ax.text(0.03, 0.70, "全样本 rho = 0.917", transform=ax.transAxes,
        fontsize=5.5, va="top", color="#333333")
ax.set_title("层内相关（分层后仍成立）", pad=5)
S.ygrid(ax); S.clean(ax); S.panel(ax, "f", x=-0.30, y=1.03)

S.save(fig, "Fig6_perturbation")