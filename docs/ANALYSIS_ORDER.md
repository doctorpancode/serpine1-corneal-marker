# Analysis order

The study ran as nine daily batches (D1–D9). Each directory keeps the scripts that produced that
day's tables; scripts are numbered in execution order, and `b`/`c` suffixes are fix-ups kept for
provenance (the final state of each analysis is the highest suffix).

## 0. Inputs

Place the downloaded matrices in `data/raw/`:

```
data/raw/
  GSE275678/            mouse scRNA-seq (alkali burn), matrix + barcodes + features
  GSE241715/            human bulk microarray series matrix
  GSE58291_expr.tsv     Illumina HT-12 v4, probe / symbol / 30 *.AVG_Signal columns
  GSE294213/ GSE204791/ GSE89537/ GSE326016/ GSE319836/   series matrices
```

`GSE58291_expr.tsv` keeps the original column names (`M 23.AVG_Signal`); the loader normalises them
to `M23` by removing whitespace — do not strip leading letters, or sample names collapse to bare
numbers and the group join silently fails.

## 1. Per-day chain

| Day | Entry script | Consumes | Produces |
|---|---|---|---|
| D2 | `D2/01_build_stemness_geneset.R` → `02_stemness_wgcna.R` → `03_export_module_genes.R` → `04_module_enrichment.R` → `05_crosscohort_stemness.R` | GSE275678 pseudobulk | `stemness_geneset_source_log.tsv`, `wgcna_module_trait_corr.tsv`, `D2_module_enrichment_*.tsv` |
| D3 | `D3/01_candidate_pool.py` → `02_pool_restrict_limma.R` → `03_venn_plot.R` | WGCNA modules + MSigDB sets + GSE241715 | `D3_limma_deg.tsv`, `D3_pool_tiering.tsv`, `D3_venn_counts.tsv` |
| D4 | `D4/01_string_ppi.py` → `03_ml_screen.R` → `04_converge.R` → `09c_crosscohort_final.R` → `11_scorecard.R` | 167-gene pool | `D4_final_genes.tsv` (10 core genes), `D4_ml_vote_matrix.tsv` |
| D5 | `D5/01_bulk_portrait.R` → `02_external_validation.R` → `03_logistic_nomogram.R` → `04_immune_infiltration.R` → `05_enrichment.R` | GSE241715 | `D5_model_ROCAUC.tsv`, `D5_main_model_OR.tsv`, `D5_vif.tsv`, `D5_enrich_*.tsv` |
| D6 | `D6/01_build_seurat.R` → `03c_annotate_fix.R` → `04_aucell_localize.R` → `06_cell_communication.R` | GSE275678 single-cell | `D6_umap_meta.tsv` (also in `data/metadata/`), `D6_aucell_summary.tsv`, `D6_celltype_composition_wide.tsv` |
| D7 | `D7/01_deconvolution.R` → `02c_sctenifold_official.R` → `02d_sctenifold_two_ko.R` → `03b_enrichment_official.R` → `04b_main_figures_official.R` | 494-gene network | `D7_deconvolution_by_group.tsv`, `D7_knk_official_twoKO.tsv`, `D7_coregene_vs_celltype_cor.tsv` |
| D8 | `D8/01_build_network.R` → `02_transcriptome_wide_scan.R` → `03_third_ko.R` → `04_spectrum_figure.R` | 947-gene network | `D8_third_ko.tsv`, `D8_knockdown_spectrum.tsv`, `D8_epithelial_program_score.tsv`, `D8_epi_selectivity_ratio.tsv` |
| D9 | `D9/01_inflammation_validation.R` → `02_selectivity_stats_and_centrality_fig.R` → `03_regularized_knockout.R` → `04_multiko.R` → `05_deconvolution.R` → `06_real_cellchat.R` → `07_crossdataset.R` | D8 outputs + external cohorts | `D9_core_genes_GSE58291.tsv`, `D9_core_genes_residual_rank.tsv`, `D9_multiko_summary.tsv`, `D9_regularized_ko_scores.tsv`, `D9_cellchat_*` |

## 2. Manuscript layer

```
code/manuscript/rev12_analysis.py     revision statistics -> rev12_stats.json
code/manuscript/build_draft12.py      assembles 草稿1.2_带图.docx (needs the template docx)
code/manuscript/dump_draft12.py       flattens the docx to text for number auditing
code/figures/fig1..fig7, s1..s8, table1   figures; read results/ through figstyle.py
```

`figstyle.py` resolves paths relative to the repository root (override with `VKO_BASE`), and writes
PNGs/PDFs into `figures_out/` (override with `VKO_OUT`).

## 3. R environment

```r
# core
install.packages(c("WGCNA", "limma", "GSVA", "preprocessCore", "clusterProfiler",
                   "org.Mm.eg.db", "org.Hs.eg.db", "e1071", "randomForest", "glmnet"))
# single-cell
remotes::install_github("satijalab/seurat", ref = "v5.1.0")
remotes::install_github("jinworks/CellChat")
# virtual knockout
remotes::install_github("dcjones/scTenifoldKnk")
```

## 4. Determinism

`scTenifoldKnk` and the network constructions are stochastic. Figures were produced from the stored
result tables, not by re-running the networks, so the shipped numbers are the ones the manuscript
reports. Re-running D7/D8 from scratch will give slightly different edge weights; the ranking claims
in the paper are stated as point estimates for exactly this reason.