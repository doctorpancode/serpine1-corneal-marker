# ============================================================
# D5-04 免疫浸润（ssGSEA, 28 细胞类型）+ 免疫检查点 + 主干基因关联
# 细胞类型标记集：Charoentong et al. 2017 (Cell Rep) 28 immune cell types
# 输出：results/D5_immune_ssgsea.tsv, D5_immune_group_stats.tsv,
#       results/D5_immune_correlation.tsv, D5_checkpoint_expression.tsv
#       figures/D5_Fig8_immune_infiltration.pdf/.png, D5_Fig9_checkpoint.pdf/.png
# ============================================================
suppressPackageStartupMessages({
  library(Biobase); library(data.table); library(GSVA); library(ggplot2)
  library(patchwork)
})
out <- "D:/A生信分析/虚拟敲除/第5天执行结果"
DATA <- "D:/A生信分析/虚拟敲除/D1执行结果/D1_CNV_角膜新生血管_数据分析/data/processed"

# Charoentong 28 immune cell signature（每类 6–15 个代表基因）
sig <- list(
  Activated_B          = c("CD79A","CD79B","MS4A1","TNFRSF17","CD19"),
  Activated_CD4        = c("CD28","ICOS","TNFRSF9","CD40LG","IL2RA"),
  Activated_CD8        = c("CD8A","CD8B","GZMK","CD27","IFNG"),
  Activated_DC         = c("CD83","CCR7","LAMP3","CD40","BIRC3"),
  CD56bright_NK        = c("NCAM1","KIR2DL4","GZMK","CD96","SELL"),
  CD56dim_NK           = c("KIR2DL3","KIR3DL1","FCGR3A","GZMB","PRF1"),
  Eosinophil           = c("SIGLEC8","CCR3","PRG2","IL5RA","EPX"),
  Gamma_delta_T        = c("TRDC","TRGC1","TRGC2","KLRD1","CD3D"),
  Immature_B           = c("IGHM","IGHD","CD38","VPREB3","DNTT"),
  Immature_DC          = c("CD1C","FCER1A","CLEC10A","CD1E","CSF1R"),
  Macrophage           = c("CD68","CD163","MRC1","CSF1R","MSR1"),
  Mast                 = c("TPSAB1","TPSB2","CPA3","MS4A2","KIT"),
  MDSC                 = c("ARG1","S100A8","S100A9","ITGAM","CD33"),
  Memory_B             = c("CD27","AICDA","IGHG1","IGHG3","MEF2B"),
  Monocyte             = c("CD14","FCGR3A","LYZ","CSF1R","S100A12"),
  Neutrophil           = c("FCGR3B","CSF3R","CXCR2","ELANE","MPO"),
  NK                   = c("KLRD1","KLRF1","NCR1","SH2D1B","NCAM1"),
  NK_T                 = c("CD3D","CD3E","KLRB1","NCAM1","TRDC"),
  NKT                  = c("CD3D","CD3E","KLRB1","ZBTB16","NCAM1"),
  Plasmacytoid_DC      = c("IL3RA","IRF7","GZMB","CLEC4C","TCF4"),
  Regulatory_T         = c("FOXP3","IL2RA","CTLA4","IKZF2","TNFRSF18"),
  T_follicular_helper  = c("CXCR5","ICOS","BCL6","IL21","PDCD1"),
  Type1_T_helper       = c("TBX21","IFNG","CXCR3","IL12RB2","STAT4"),
  Type17_T_helper      = c("RORC","IL17A","IL23R","CCR6","IL17F"),
  Type2_T_helper       = c("GATA3","IL4","IL5","IL13","CCR4"),
  Central_memory       = c("CCR7","SELL","IL7R","TCF7","CD27"),
  Effector_memory      = c("GZMK","CD44","CD27","KLRG1","EOMES"),
  Activated_B2         = c("BANK1","FCRL3","IGHG2","POU2AF1","SPIB")
)

es <- readRDS(file.path(DATA, "GSE241715_expression.rds"))
x <- exprs(es); pd <- pData(es)
sig <- lapply(sig, function(g) intersect(g, rownames(x)))
sig <- sig[lengths(sig) >= 2]
cat("可用细胞类型:", length(sig), "\n")

gp <- ssgseaParam(x, sig, normalize = TRUE, minSize = 2, maxSize = 200)
ss <- gsva(gp, verbose = FALSE)
fwrite(data.table(celltype = rownames(ss), ss,
                  check.names = FALSE), file.path(out, "results", "D5_immune_ssgsea.tsv"), sep = "\t")

d <- data.table(t(ss))
d[, group := factor(pd$group, levels = c("Control", "Keratoconus", "Viral", "Bacterial_fungal"))]
d[, lesion := ifelse(pd$group == "Control", "Normal", "Lesion")]
stat <- rbindlist(lapply(rownames(ss), function(ct) {
  v <- as.numeric(ss[ct, ])
  p <- suppressWarnings(wilcox.test(v ~ d$lesion)$p.value)
  k <- suppressWarnings(kruskal.test(v ~ d$group)$p.value)
  data.table(celltype = ct,
             mean_normal = mean(v[d$lesion == "Normal"]),
             mean_lesion = mean(v[d$lesion == "Lesion"]),
             delta = mean(v[d$lesion == "Lesion"]) - mean(v[d$lesion == "Normal"]),
             p_wilcox = p, p_kruskal = k)
}))
stat[, adjP := p.adjust(p_wilcox, "BH")]
setorder(stat, p_wilcox)
fwrite(stat, file.path(out, "results", "D5_immune_group_stats.tsv"), sep = "\t")
print(stat[, .(celltype, delta = round(delta, 3), p = signif(p_wilcox, 3), adjP = signif(adjP, 3))])

# ---- 主干基因与免疫细胞相关性 ----
corex <- t(x[c("TMSB4X", "SERPINE1"), , drop = FALSE])
corl <- rbindlist(lapply(rownames(ss), function(ct) {
  rbindlist(lapply(c("TMSB4X", "SERPINE1"), function(g) {
    ct2 <- cor.test(corex[, g], ss[ct, ], method = "spearman")
    data.table(gene = g, celltype = ct, rho = as.numeric(ct2$estimate), p = ct2$p.value)
  }))
}))
corl[, adjP := p.adjust(p, "BH")]
setorder(corl, p)
fwrite(corl, file.path(out, "results", "D5_immune_correlation.tsv"), sep = "\t")
cat("\n=== SERPINE1 相关 top ===\n"); print(corl[gene == "SERPINE1"][1:8])
cat("\n=== TMSB4X 相关 top ===\n"); print(corl[gene == "TMSB4X"][1:8])

# ---- 免疫检查点 ----
ck <- c("PDCD1","CD274","PDCD1LG2","CTLA4","HAVCR2","TIGIT","LAG3","VSIR","BTLA","IDO1","CD80","CD86")
ck <- intersect(ck, rownames(x))
ckd <- data.table(t(x[ck, , drop = FALSE]))
ckd[, group := factor(pd$group, levels = c("Control", "Keratoconus", "Viral", "Bacterial_fungal"))]
ckd[, lesion := ifelse(pd$group == "Control", "Normal", "Lesion")]
ckstat <- rbindlist(lapply(ck, function(g) {
  v <- ckd[[g]]
  data.table(checkpoint = g,
             mean_normal = mean(v[ckd$lesion == "Normal"]),
             mean_lesion = mean(v[ckd$lesion == "Lesion"]),
             p_wilcox = suppressWarnings(wilcox.test(v ~ ckd$lesion)$p.value))
}))
ckstat[, adjP := p.adjust(p_wilcox, "BH")]
setorder(ckstat, p_wilcox)
fwrite(ckstat, file.path(out, "results", "D5_checkpoint_expression.tsv"), sep = "\t")
print(ckstat)

# ---- 图8：免疫浸润热图（z-score of ssGSEA）+ 差异条形 ----
suppressPackageStartupMessages({library(ComplexHeatmap); library(circlize)})
m <- t(scale(t(ss)))
gcol <- c(Control = "#4c72b0", Keratoconus = "#55a868", Viral = "#dd8452", Bacterial_fungal = "#d62728")
ha <- HeatmapAnnotation(Group = as.character(pd$group), col = list(Group = gcol),
                        annotation_name_side = "left")
ht <- Heatmap(m, name = "ssGSEA\nz-score",
              col = colorRamp2(c(-2, 0, 2), c("#2166ac", "white", "#b2182b")),
              top_annotation = ha, cluster_rows = TRUE, cluster_columns = TRUE,
              row_names_side = "left", row_names_gp = gpar(fontsize = 7),
              column_names_gp = gpar(fontsize = 7),
              column_title = "Immune infiltration (28 cell types, ssGSEA; GSE241715)",
              column_title_gp = gpar(fontsize = 9))
pdf(file.path(out, "figures", "D5_Fig8_immune_infiltration.pdf"), width = 7.6, height = 6.2)
draw(ht); dev.off()
png(file.path(out, "figures", "D5_Fig8_immune_infiltration.png"), width = 7.6, height = 6.2,
    units = "in", res = 300); draw(ht); dev.off()

top <- stat[order(-abs(delta))][1:12]
pl <- melt(d[, c("group", top$celltype), with = FALSE], id.vars = "group",
           variable.name = "celltype", value.name = "score")
pI <- ggplot(pl, aes(group, score, fill = group)) +
  geom_boxplot(width = 0.62, outlier.size = 0.7) +
  facet_wrap(~ celltype, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = gcol) +
  labs(title = "Top-12 immune cell types by absolute delta (Lesion vs Normal)",
       x = NULL, y = "ssGSEA score") +
  theme_bw(base_size = 9) + theme(legend.position = "none",
                                  axis.text.x = element_text(angle = 40, hjust = 1))
ggsave(file.path(out, "figures", "D5_Fig8B_immune_boxplot.pdf"), pI, width = 9.6, height = 5.4)
ggsave(file.path(out, "figures", "D5_Fig8B_immune_boxplot.png"), pI, width = 9.6, height = 5.4, dpi = 300)

# ---- 图9：检查点 + 主干基因相关 ----
ckl <- melt(ckd[, c("group", ck), with = FALSE], id.vars = "group",
            variable.name = "checkpoint", value.name = "expr")
pK <- ggplot(ckl, aes(group, expr, fill = group)) +
  geom_boxplot(width = 0.6, outlier.size = 0.6) +
  facet_wrap(~ checkpoint, scales = "free_y", ncol = 6) +
  scale_fill_manual(values = gcol) +
  labs(title = "Immune checkpoint expression across diagnostic groups (GSE241715)",
       x = NULL, y = "log2(CPM + 1)") +
  theme_bw(base_size = 8.5) + theme(legend.position = "none",
                                    axis.text.x = element_text(angle = 45, hjust = 1))
pK2 <- ggplot(corl, aes(reorder(celltype, rho), rho, fill = gene)) +
  geom_col(position = position_dodge(0.7), width = 0.68) +
  geom_hline(yintercept = 0, colour = "grey40") +
  coord_flip() +
  scale_fill_manual(values = c(TMSB4X = "#4c72b0", SERPINE1 = "#d62728")) +
  labs(title = "Spearman correlation: core genes vs immune infiltration",
       x = NULL, y = "rho", fill = NULL) +
  theme_bw(base_size = 8.5) + theme(legend.position = "top")
ggsave(file.path(out, "figures", "D5_Fig9_checkpoint.pdf"),
       wrap_plots(pK, pK2, ncol = 2, widths = c(1.45, 1)), width = 13, height = 5.6)
ggsave(file.path(out, "figures", "D5_Fig9_checkpoint.png"),
       wrap_plots(pK, pK2, ncol = 2, widths = c(1.45, 1)), width = 13, height = 5.6, dpi = 300)
saveRDS(ss, file.path(out, "data", "D5_immune_ssgsea.rds"))
cat("\nD5-04 done\n")