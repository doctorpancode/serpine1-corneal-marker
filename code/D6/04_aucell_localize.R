# ============================================================
# D6-04 AUCell 定位两个主干基因（Serpine1 / Tmsb4x）
# 同时给出所有候选基因在各细胞类型的平均表达与检出率
# ============================================================
suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(ggplot2); library(patchwork); library(AUCell)
})
set.seed(20261109)
out <- "D:/A生信分析/虚拟敲除/第6天执行结果"; fig <- file.path(out, "figures")

o <- readRDS(file.path(out, "data", "D6_seurat_GSE275678_final.rds"))
cand <- c("Serpine1","Tmsb4x","P4hb","Vdac1","Sod1","Prdx2","Gja1","S100a13","Capg","Xiap")
cand <- intersect(cand, rownames(o))
cat("candidates present:", paste(cand, collapse = ", "), "\n")

# ---- 1) 平均表达 + 检出率（按 celltype）----
expr <- GetAssayData(o, layer = "data")
ct <- o$celltype
res <- rbindlist(lapply(cand, function(g) {
  v <- expr[g, ]
  dt <- data.table(gene = g, celltype = ct, expr = as.numeric(v))
  dt[, .(mean_expr = mean(expr), pct_expr = 100 * mean(expr > 0), n = .N), by = celltype]
}))
fwrite(res, file.path(out, "results", "D6_candidate_expression_by_celltype.tsv"), sep = "\t")

# 按 sample 再分
res2 <- rbindlist(lapply(cand, function(g) {
  v <- as.numeric(expr[g, ])
  dt <- data.table(gene = g, celltype = ct, sample = o$sample, expr = v)
  dt[, .(mean_expr = mean(expr), pct_expr = 100 * mean(expr > 0), n = .N), by = .(celltype, sample)]
}))
fwrite(res2, file.path(out, "results", "D6_candidate_expression_celltype_sample.tsv"), sep = "\t")

# ---- 2) 点图：两主干在细胞类型 × 样本 ----
o$grp <- paste0(o$celltype, "|", o$sample)
dp <- DotPlot(o, features = c("Serpine1","Tmsb4x"), group.by = "grp") +
  RotatedAxis() + ggtitle("Serpine1 / Tmsb4x across cell types and conditions")
ggsave(file.path(fig, "D6_Fig5_core_gene_dotplot.png"), dp, width = 8, height = 8, dpi = 300)
ggsave(file.path(fig, "D6_Fig5_core_gene_dotplot.pdf"), dp, width = 8, height = 8)

# ---- 3) AUCell（每基因单基因集，评估表达活性）----
rankings <- AUCell_buildRankings(as.matrix(expr), nCores = 4, plotStats = FALSE, verbose = FALSE)
auc_list <- lapply(cand, function(g) {
  s <- list(g); names(s) <- g
  a <- AUCell_calcAUC(s, rankings, aucMaxRank = ceiling(0.05 * nrow(rankings)), verbose = FALSE)
  data.table(cell = colnames(a), gene = g, auc = as.numeric(getAUC(a)))
})
auc <- rbindlist(auc_list)
auc[, celltype := o$celltype[cell]]; auc[, sample := o$sample[cell]]
fwrite(auc, file.path(out, "results", "D6_aucell_scores_long.tsv"), sep = "\t")
agg <- auc[, .(mean_auc = mean(auc), med_auc = median(auc), n = .N), by = .(gene, celltype, sample)]
fwrite(agg, file.path(out, "results", "D6_aucell_summary.tsv"), sep = "\t")
print(dcast(agg[gene == "Serpine1"], celltype ~ sample, value.var = "mean_auc"))
print(dcast(agg[gene == "Tmsb4x"], celltype ~ sample, value.var = "mean_auc"))

p_auc <- ggplot(agg, aes(x = reorder(celltype, mean_auc), y = mean_auc, fill = sample)) +
  geom_col(position = "dodge") + coord_flip() + facet_wrap(~ gene, scales = "free_x") +
  scale_fill_manual(values = c(uninjury = "#4C72B0", NaOH = "#C44E52")) +
  labs(x = "", y = "Mean AUCell score", title = "Regulon activity of core genes by cell type") +
  theme_bw(base_size = 11)
ggsave(file.path(fig, "D6_Fig6_aucell_barplot.png"), p_auc, width = 13, height = 6, dpi = 300)
ggsave(file.path(fig, "D6_Fig6_aucell_barplot.pdf"), p_auc, width = 13, height = 6)

# ---- 4) UMAP 上叠加表达 ----
pA <- FeaturePlot(o, features = "Serpine1", reduction = "umap", pt.size = 0.3) +
  scale_colour_gradientn(colours = c("grey90","#FEE08B","#F46D43","#A50026")) + ggtitle("Serpine1")
pB <- FeaturePlot(o, features = "Tmsb4x", reduction = "umap", pt.size = 0.3) +
  scale_colour_gradientn(colours = c("grey90","#D9F0D3","#66BD63","#005A32")) + ggtitle("Tmsb4x")
ggsave(file.path(fig, "D6_Fig7_core_gene_umap.png"), wrap_plots(pA, pB, ncol = 2), width = 13, height = 5.5, dpi = 300)
ggsave(file.path(fig, "D6_Fig7_core_gene_umap.pdf"), wrap_plots(pA, pB, ncol = 2), width = 13, height = 5.5)

cat("\nD6-04 done\n")