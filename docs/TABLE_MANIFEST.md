# Number → file manifest

Every quantitative claim in the manuscript maps to a file in `results/`. Files are grouped by the
analysis day that produced them (`results/D2` … `results/D9`).

## Table 1 — discriminative performance (discovery vs external cohort)

| Reported | Value | File | Column |
|---|---|---|---|
| SERPINE1 AUC, discovery | 0.949 (0.837–1.000) | `results/D5/D5_model_ROCAUC.tsv` | `AUC`, `CI_low`, `CI_high` for `M2_SERPINE1` |
| Hosmer–Lemeshow P | 0.805 | `results/D5/D5_model_ROCAUC.tsv` | `HL_p` |
| SERPINE1 AUC, LOOCV | 0.878 | `code/manuscript/rev12_analysis.py` (recomputed, 2 000 stratified bootstrap) | — |
| Perfect/near-perfect AUCs | 1.000 (P4HB, S100A13, CAPG, XIAP); 0.000 (TMSB4X, GJA1, SOD1); 0.010 (PRDX2, VDAC1) | `code/manuscript/rev12_analysis.py` | — |
| SERPINE1 AUC, external | 0.988 (0.929–1.000) | `code/manuscript/rev12_analysis.py` | — |
| SERPINE1 AUC, infections combined | 0.956 (0.856–1.000) | `code/manuscript/rev12_analysis.py` | — |
| TMSB4X AUC, external | 0.488 (0.214–0.762) | `code/manuscript/rev12_analysis.py` | — |
| Univariate OR | not estimable (complete separation) | `results/D5/D5_main_model_OR.tsv` | `beta`, `OR`, `CI_low`, `CI_high`, `p` |
| VIF | all > 5; GJA1 70.5 … PRDX2 9.0 | `results/D5/D5_vif.tsv` | `VIF` |

## Table S1 — 18 epithelial genes under TMSB4X / SERPINE1 knockout

| Reported | File |
|---|---|
| Signed Δ per gene (18 rows) | `data/metadata/D8_epi18_delta_subset.tsv` (subset of `D8_knockdown_spectrum.tsv`, 54 MB, in the Zenodo archive only) |
| mean\|Δ\| across the 18 genes | `results/D8/D8_epithelial_program_score.tsv` |

## Table S2 — key numbers

| Reported | Value | File |
|---|---|---|
| DEGs, FDR < 0.05 | 1 367 | `results/D3/D3_limma_deg.tsv` (row count with `adjP < 0.05`) |
| DEGs, FDR < 0.05 and \|log2FC\| > 0.5 | 1 350 (743 ↑ / 607 ↓) | same file |
| SERPINE1 group means, discovery | 1.648 / 3.330 / 5.226 / 8.904 | `code/manuscript/rev12_analysis.py` on `data/metadata/GSE241715_sample_sheet.tsv` + GEO matrix |
| Normal → bacterial/fungal fold change | 5.40 | same |
| SERPINE1 means, external | 442.9 / 1 733.3 / 4 777.6 (log2 8.48 / 10.50 / 11.54) | `results/D9/D9_core_genes_GSE58291.tsv` |
| BK vs control | log2FC +3.06, adjP 2.5×10⁻⁷, ρ = 0.828 | `results/D9/D9_core_genes_GSE58291.tsv` |
| Epithelial mean\|Δ\| TMSB4X / P4HB / SERPINE1 | 0.271213 / 0.008837 / 0.000587 (462×) | `results/D9/D9_multiko_summary.tsv`, `results/D8/D8_knockdown_spectrum.tsv` |
| Largest epithelial effect of SERPINE1 | 1.587×10⁻³ (Krt5) | `results/D8/D8_knockdown_spectrum.tsv` |
| Residual rank after out-degree detrending | TMSB4X 1 (0.246173), SERPINE1 750 (−0.002386), next 0.018806 | `results/D9/D9_core_genes_residual_rank.tsv` |
| Out-degree–epithelial correlation | ρ 0.9169; within-tercile 0.9498 / 0.4909 / 0.6055 | `results/D8/D8_ko_impact_rank.tsv`, `results/D9/D9_regularized_ko_scores.tsv` |
| Row-L1 renormalised rank shift | Serpine1 648→764, Prdx2 50→14, Vdac1 175→276, Tmsb4x 1→1 | `results/D9/D9_regularized_ko_scores.tsv` |
| Multi-knockout ratios | 1.995 / 1.893 / 1.874; triple 2.805; TMSB4X+P4HB mean\|Δ\| 0.262472 < single 0.271213 | `results/D9/D9_multiko_summary.tsv` |
| Non-target significant genes (5 knockouts) | 0 (minimum adjP 0.438–0.992) | `results/D8/D8_third_ko_summary.tsv`, `results/D8/D8_third_ko.tsv` |
| Epithelial selectivity test | median 1.348 vs 1.150; Wilcoxon P 0.6755; permutation P 0.4130; δ −0.0796 | `results/D8/D8_epi_selectivity_ratio.tsv` |
| Network size | 494 × 2 006 (169 562 edges, 47.7 %); 947 × 3 587 (600 604 edges, 66.97 %, CP K=3 56.4 %) | day 7 / day 8 logs (`code/D7/01_deconvolution.R`, `code/D8/01_build_network.R`) |
| Cross-model direction concordance | 2 models 9/10; 3 models 3/10 (Vdac1, Capg, Serpine1) | `results/D9/D9_core_gene_direction_concordance_3ds.tsv` |
| Cell-type composition | neutrophil 0.029 %→11.335 %, macrophage 0.245 %→4.330 %, basal epithelium 43.8 %→21.8 %, fibroblast 1 cell injured side | `results/D6/D6_celltype_composition_wide.tsv`, `data/metadata/D6_umap_meta.tsv` |
| SERPINE1 AUCell activity | endothelium 0.00374→0.22339 (59.7×) | `results/D6/D6_aucell_summary.tsv` |
| SERPINE1–cell-type correlation | macrophage ρ 0.821, epithelial wing ρ −0.764 | `results/D7/D7_coregene_vs_celltype_cor.tsv` |
| Deconvolution weights | 4 groups × 11 cell types | `results/D7/D7_deconvolution_by_group.tsv` |
| CellChat | 65 interactions, SPP1 0→0.0180, CCL 2 854× | `results/D9/D9_cellchat_compare_pathways.tsv`, `results/D9/D9_cellchat_lr_to_vessel_real.tsv` |

## Notes on provenance

- `rev12_analysis.py` produced the numbers added during revision (discovery-cohort per-sample values,
  LOOCV, external AUCs). It reads only public GEO matrices plus the sample sheets shipped in
  `data/metadata/`; its JSON output is `rev12_stats.json` (regenerate on demand).
- The four matrices larger than 2 MB are not in the GitHub tree. `D8_epi18_delta_subset.tsv` is
  shipped so that Table S1 can be reproduced without them.
- `D3_limma_deg.tsv` has no contrast column (it holds a single lesion-vs-control fit), which is why
  the DEG counts are obtained by row counting rather than by grouping.