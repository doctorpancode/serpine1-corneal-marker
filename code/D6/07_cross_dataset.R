# ============================================================
# D6-07 跨物种/跨数据集一致性（GSE294213 缝线 CNV 验证）
# 用同一套 marker 注释 + 检查 Serpine1/Tmsb4x 在两队列中的方向
# ============================================================
suppressPackageStartupMessages({ library(Seurat); library(data.table); library(ggplot2) })
set.seed(20261109)
out <- "D:/A生信分析/虚拟敲除/第6天执行结果"; fig <- file.path(out, "figures")

o <- readRDS(file.path(out, "data", "D6_seurat_GSE294213_qc.rds"))
cat("GSE294213 cells:", ncol(o), " genes:", nrow(o), "\n")

o <- NormalizeData(o, verbose = FALSE)
o <- FindVariableFeatures(o, nfeatures = 2000, verbose = FALSE)
o <- ScaleData(o, verbose = FALSE)
o <- RunPCA(o, npcs = 30, verbose = FALSE)
o <- FindNeighbors(o, dims = 1:30, verbose = FALSE)
o <- FindClusters(o, resolution = 0.5, verbose = FALSE)
o <- RunUMAP(o, dims = 1:30, verbose = FALSE)
o <- JoinLayers(o)
cat("clusters:", length(levels(o$seurat_clusters)), "\n")

mk <- as.data.table(FindAllMarkers(o, only.pos = TRUE, min.pct = 0.25,
                                   logfc.threshold = 0.5, verbose = FALSE))
fwrite(mk, file.path(out, "results", "D6_GSE294213_top_markers.tsv"), sep = "\t")
top3 <- mk[order(cluster, -avg_log2FC)][, head(.SD, 3), by = cluster]
print(top3[, .(cluster, gene, avg_log2FC = round(avg_log2FC, 2))])

# marker 打分注释
marker_sets <- list(
  Epithelium   = c("Krt3","Krt12","Krt5","Krt14","Krt15","Trp63","Krt19","Krt13"),
  Endothelium  = c("Pecam1","Vwf","Cdh5","Kdr","Emcn","Esam"),
  Lymphatic    = c("Lyve1","Prox1","Pdpn","Ccl21a","Flt4"),
  Fibroblast   = c("Col1a1","Col1a2","Col3a1","Lum","Dcn","Pdgfra"),
  Macrophage   = c("Cd68","Adgre1","Csf1r","Lgals3","Itgam"),
  Neutrophil   = c("S100a8","S100a9","Ly6g6c","Acod1","Cxcr2"),
  Tcell        = c("Cd3e","Cd3d","Trac","Lck","Cd2"),
  Pericyte_SMC = c("Rgs5","Kcnj8","Abcc9","Myh11","Des"),
  Proliferating= c("Mki67","Top2a","Cenpf","Nusap1","Mcm3","Hist1h3c")
)
for (nm in names(marker_sets)) {
  g <- intersect(marker_sets[[nm]], rownames(o))
  if (length(g) >= 2) o <- AddModuleScore(o, features = list(g), name = nm, seed = 20261109)
}
sc_cols <- grep("1$", grep("^[A-Z]", names(o@meta.data), value = TRUE), value = TRUE)
md <- as.data.table(o@meta.data); md[, cluster := as.character(seurat_clusters)]
sc <- md[, c(list(n = .N), lapply(.SD, mean)), by = cluster, .SDcols = sc_cols]
fwrite(sc, file.path(out, "results", "D6_GSE294213_cluster_scores.tsv"), sep = "\t")
print(sc[, 1:5])

# 注意：sc 的行序是 by=cluster 排序后的顺序，必须用它自己的 cluster 列做映射，
# 且 sc_cols 是 meta.data 里的列名（带 1 后缀，如 Epithelium1）
sc_cols <- grep("^[A-Za-z_]+1$", names(sc), value = TRUE)
sc_cols <- setdiff(sc_cols, "n")
cat("score columns:", paste(sc_cols, collapse = ", "), "\n")
zsc <- scale(as.matrix(sc[, ..sc_cols]))
best <- sc_cols[apply(zsc, 1, which.max)]
ann <- setNames(sub("1$", "", best), as.character(sc$cluster))
o$celltype <- unname(ann[as.character(o$seurat_clusters)])
o$celltype[is.na(o$celltype)] <- "Unknown"
cat("\ncelltype assignment:\n"); print(table(o$celltype))
cat("all assigned:", all(!is.na(o$celltype)), "\n")
saveRDS(o, file.path(out, "data", "D6_seurat_GSE294213_final.rds"))

# ---- 两队列归一化表达对比（同基因同条件）----
o1 <- readRDS(file.path(out, "data", "D6_seurat_GSE275678_stem.rds"))
genes <- c("Serpine1","Tmsb4x","P4hb","Vdac1","Sod1","Prdx2","Gja1","S100a13","Capg","Xiap")
genes <- intersect(genes, intersect(rownames(o1), rownames(o)))
cmp <- rbindlist(lapply(genes, function(g) {
  e1 <- GetAssayData(o1, layer = "data")
  v1 <- as.numeric(e1[g, ])
  d1 <- data.table(dataset = "GSE275678", gene = g,
                   cond = ifelse(o1$sample == "NaOH", "Lesion", "Control"),
                   celltype = o1$celltype, expr = v1)
  e2 <- GetAssayData(o, layer = "data")
  v2 <- as.numeric(e2[g, ])
  d2 <- data.table(dataset = "GSE294213", gene = g,
                   cond = ifelse(o2s <- (o$sample == "nosut"), "Control", "Lesion"),
                   celltype = o$celltype, expr = v2)
  rbind(d1, d2)
}))
agg <- cmp[, .(mean_expr = mean(expr), pct = 100 * mean(expr > 0), n = .N),
           by = .(dataset, gene, cond, celltype)]
fwrite(agg, file.path(out, "results", "D6_cross_dataset_celltype_expr.tsv"), sep = "\t")

# 每个数据集的整体方向（所有细胞合并）
agg2 <- cmp[, .(mean_expr = mean(expr), n = .N), by = .(dataset, gene, cond)]
w <- dcast(agg2, dataset + gene ~ cond, value.var = "mean_expr")
w[, delta := Lesion - Control]
w[, direction := ifelse(delta > 0, "up", "down")]
fwrite(w, file.path(out, "results", "D6_cross_dataset_direction.tsv"), sep = "\t")
print(w)

p <- ggplot(w, aes(x = gene, y = delta, fill = direction)) +
  geom_col() + facet_wrap(~ dataset, scales = "free_y") + coord_flip() +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_fill_manual(values = c(up = "#C44E52", down = "#4C72B0")) +
  labs(x = "", y = "Delta mean expression (Lesion - Control)",
       title = "Cross-dataset consistency of candidate genes (scRNA)") +
  theme_bw(base_size = 12)
ggsave(file.path(fig, "D6_Fig11_cross_dataset.png"), p, width = 10, height = 5.5, dpi = 300)
ggsave(file.path(fig, "D6_Fig11_cross_dataset.pdf"), p, width = 10, height = 5.5)

cat("\nD6-07 done\n")