# ============================================================
# D6-02 单细胞标准流程 + 细胞注释（GSE275678 主分析）
# PCA -> UMAP -> cluster -> marker 注释
# ============================================================
suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(ggplot2); library(patchwork)
})
set.seed(20261109)
out <- "D:/A生信分析/虚拟敲除/第6天执行结果"
fig <- file.path(out, "figures"); dir.create(fig, recursive = TRUE, showWarnings = FALSE)

o <- readRDS(file.path(out, "data", "D6_seurat_GSE275678_qc.rds"))
o <- NormalizeData(o, verbose = FALSE)
o <- FindVariableFeatures(o, nfeatures = 2000, verbose = FALSE)
o <- ScaleData(o, verbose = FALSE)
o <- RunPCA(o, npcs = 30, verbose = FALSE)
o <- FindNeighbors(o, dims = 1:30, verbose = FALSE)
o <- FindClusters(o, resolution = 0.6, verbose = FALSE)
o <- RunUMAP(o, dims = 1:30, verbose = FALSE)

cat("clusters:", length(levels(o$seurat_clusters)), "\n")
cat("cells per cluster:\n"); print(table(o$seurat_clusters))

# ---- marker 打分（自实现：AddModuleScore 等价）----
marker_sets <- list(
  Cornea_epithelial = c("Krt3","Krt12","Krt5","Krt14","Krt15","Trp63","Krt19","Krt13"),
  Endothelial       = c("Pecam1","Vwf","Cdh5","Kdr","Emcn","Esam","Cd34"),
  Lymphatic         = c("Lyve1","Prox1","Pdpn","Ccl21a","Flt4"),
  Fibroblast        = c("Col1a1","Col1a2","Col3a1","Lum","Dcn","Pdgfra","Col5a1"),
  Myofibroblast     = c("Acta2","Tagln","Cnn1","Myl9","Transgelin"),
  Macrophage        = c("Cd68","Adgre1","Csf1r","Lgals3","Itgam","Cd14"),
  Tcell             = c("Cd3e","Cd3d","Cd3g","Trac","Lck","Cd2"),
  Bcell             = c("Cd79a","Cd79b","Ms4a1","Ighm","Cd19"),
  Neutrophil        = c("S100a8","S100a9","Ly6g","Retnlg","Mpo","Elane"),
  Pericyte_SMC      = c("Rgs5","Kcnj8","Abcc9","Myh11","Des","Notch3"),
  Melanocyte        = c("Pmel","Mlana","Tyr","Dct","Pax3"),
  Schwann_Neural    = c("Plp1","Mpz","Sox10","S100b","Ngfr")
)
for (nm in names(marker_sets)) {
  g <- intersect(marker_sets[[nm]], rownames(o))
  o <- AddModuleScore(o, features = list(g), name = nm, seed = 20261109)
  cat(nm, ":", length(g), "markers used\n")
}
score_cols <- paste0(names(marker_sets), "1")
md <- as.data.table(o@meta.data)
md[, cluster := as.character(seurat_clusters)]
sc <- md[, c(list(n = .N), lapply(.SD, mean)), by = cluster, .SDcols = score_cols]
fwrite(sc, file.path(out, "results", "D6_cluster_marker_scores.tsv"), sep = "\t")
print(sc[order(-n)][1:12])

saveRDS(o, file.path(out, "data", "D6_seurat_GSE275678_annot_input.rds"))

# ---- Fig: UMAP by cluster / by sample ----
p1 <- DimPlot(o, reduction = "umap", label = TRUE, repel = TRUE) + ggtitle("Clusters (GSE275678)") + NoLegend()
p2 <- DimPlot(o, reduction = "umap", group.by = "sample", cols = c("#4C72B0","#C44E52")) + ggtitle("Sample origin")
ggsave(file.path(fig, "D6_Fig1_UMAP_cluster_sample.png"), wrap_plots(p1, p2, ncol = 2), width = 13, height = 5.5, dpi = 300)
ggsave(file.path(fig, "D6_Fig1_UMAP_cluster_sample.pdf"), wrap_plots(p1, p2, ncol = 2), width = 13, height = 5.5)

# ---- Fig: marker dotplot ----
dp_genes <- unique(unlist(marker_sets))
dp_genes <- intersect(dp_genes, rownames(o))
p3 <- DotPlot(o, features = dp_genes, group.by = "seurat_clusters") + RotatedAxis() +
  theme(axis.text.x = element_text(size = 6))
ggsave(file.path(fig, "D6_Fig2_marker_dotplot.png"), p3, width = 16, height = 6, dpi = 300)
ggsave(file.path(fig, "D6_Fig2_marker_dotplot.pdf"), p3, width = 16, height = 6)

cat("\nD6-02 done\n")