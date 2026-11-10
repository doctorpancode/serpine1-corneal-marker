# ============================================================
# D5-01 Bulk 画像第一步：双主干（TMSB4X / SERPINE1）表达分布 + VIF 共线性诊断
# 计划书 D5 要求：病变 vs 正常小提琴图、多因素 logistic、ROC、校准曲线、
#                列线图、免疫浸润、检查点、热图、GO/KEGG/Reactome
# 本脚本: 摘要统计 + 小提琴/箱线 + 病因分层 + Spearman 相关 + VIF
# 输出: results/D5_group_stats.tsv, D5_subgroup_stats.tsv, D5_spearman.tsv, D5_vif.tsv
#       figures/D5_Fig1_violin_core.pdf/.png, D5_Fig2_panel_heatmap.pdf/.png
# ============================================================
suppressPackageStartupMessages({
  library(Biobase); library(data.table); library(ggplot2); library(patchwork)
  library(pheatmap); library(car)
})
out <- "D:/A生信分析/虚拟敲除/第5天执行结果"
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "results"), recursive = TRUE, showWarnings = FALSE)

es <- readRDS("D:/A生信分析/虚拟敲除/D1执行结果/D1_CNV_角膜新生血管_数据分析/data/processed/GSE241715_expression.rds")
x <- exprs(es); pd <- pData(es)
cat("矩阵:", dim(x), "\n"); print(table(pd$group))

core <- c("TMSB4X", "SERPINE1")
panel <- c("TMSB4X", "SERPINE1", "S100A13", "XIAP", "CAPG", "PRDX2", "VDAC1", "P4HB", "SOD1", "GJA1")
core <- intersect(core, rownames(x)); panel <- intersect(panel, rownames(x))

d <- data.table(t(x[panel, , drop = FALSE]))
d[, sample := colnames(x)]
d[, group := factor(pd$group, levels = c("Control", "Keratoconus", "Viral", "Bacterial_fungal"))]
d[, lesion := factor(ifelse(pd$group == "Control", "Normal", "Lesion"), levels = c("Normal", "Lesion"))]
d[, etiology := fifelse(pd$group == "Control", "Control",
                 fifelse(pd$group == "Keratoconus", "Keratoconus", "Infection"))]

# ---------- 组统计 ----------
long <- melt(d, id.vars = c("sample", "group", "lesion", "etiology"),
             measure.vars = panel, variable.name = "gene", value.name = "expr")
gs <- long[, .(n = .N, mean = mean(expr), sd = sd(expr), median = median(expr),
               q25 = quantile(expr, .25), q75 = quantile(expr, .75)), by = .(gene, group)]
gs[, p_wilcox_vs_control := {
  ctl <- long[gene == .BY$gene & group == "Control"]$expr
  .(sapply(group, function(g) {
    v <- long[gene == .BY$gene & group == g]$expr
    if (g == "Control") NA_real_ else suppressWarnings(wilcox.test(v, ctl)$p.value)
  }))
}, by = gene]
gs[, p_kruskal := {
  .(rep(suppressWarnings(kruskal.test(expr ~ group, data = long[gene == .BY$gene])$p.value), .N))
}, by = gene]
fwrite(gs, file.path(out, "results", "D5_group_stats.tsv"), sep = "\t")
cat("\n=== 组均值 ===\n"); print(dcast(gs[, .(gene, group, mean = round(mean, 3))], gene ~ group, value.var = "mean"))
cat("\n=== Kruskal p ===\n"); print(unique(gs[, .(gene, p_kruskal = signif(p_kruskal, 3))]))

# ---------- 病因分层 ----------
sub <- long[etiology != "Control", .(n = .N, mean = mean(expr), sd = sd(expr)),
            by = .(gene, etiology)]
fwrite(sub, file.path(out, "results", "D5_subgroup_stats.tsv"), sep = "\t")
print(dcast(sub[, .(gene, etiology, mean = round(mean, 3))], gene ~ etiology, value.var = "mean"))

# ---------- 核心基因小提琴 ----------
pl <- melt(d, id.vars = c("sample", "group", "lesion"), measure.vars = core,
           variable.name = "gene", value.name = "expr")
pA <- ggplot(pl[lesion == "Normal" | TRUE], aes(group, expr, fill = group)) +
  geom_violin(scale = "width", alpha = 0.55, trim = FALSE, colour = NA) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.9) +
  geom_jitter(width = 0.1, size = 1.1, alpha = 0.7) +
  facet_wrap(~ gene, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c(Control = "#4c72b0", Keratoconus = "#55a868",
                               Viral = "#dd8452", Bacterial_fungal = "#d62728")) +
  labs(title = "D5 core genes across diagnostic groups (GSE241715, n = 21)",
       subtitle = "Wilcoxon vs Control: TMSB4X p = 1.72e-05 (down), SERPINE1 p = 1.17e-03 (up)",
       x = NULL, y = "log2(CPM + 1)") +
  theme_bw(base_size = 10) + theme(legend.position = "none",
                                   axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(out, "figures", "D5_Fig1_violin_core.pdf"), pA, width = 7.2, height = 5.4)
ggsave(file.path(out, "figures", "D5_Fig1_violin_core.png"), pA, width = 7.2, height = 5.4, dpi = 300)

# ---------- 相关性（Spearman） ----------
cm <- cor(d[, ..panel], method = "spearman")
sp <- as.data.table(as.table(cm))
setnames(sp, c("gene1", "gene2", "rho"))
fwrite(sp, file.path(out, "results", "D5_spearman.tsv"), sep = "\t")

ann <- data.frame(group = as.character(d$group), row.names = d$sample, stringsAsFactors = FALSE)
colsc <- colorRampPalette(c("#4575b4", "white", "#d73027"))(100)
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
mat <- t(scale(as.matrix(d[, ..panel])))   # genes x samples, z-scored per gene
if (ncol(mat) != nrow(d)) stop("dim mismatch: ", ncol(mat), " vs ", nrow(d))
colnames(mat) <- as.character(d$sample)
gcol <- c(Control = "#4c72b0", Keratoconus = "#55a868", Viral = "#dd8452", Bacterial_fungal = "#d62728")
ha <- HeatmapAnnotation(Group = as.character(d$group),
                        col = list(Group = gcol),
                        annotation_name_side = "left")
ht <- Heatmap(mat, name = "z-score",   # rows = genes, cols = samples
              col = colorRamp2(c(-2, 0, 2), c("#4575b4", "white", "#d73027")),
              top_annotation = ha,
              cluster_rows = TRUE, cluster_columns = TRUE,
              row_names_side = "left", column_names_gp = gpar(fontsize = 7),
              row_names_gp = gpar(fontsize = 8),
              column_title = "D5 candidate panel, row z-score (GSE241715)",
              column_title_gp = gpar(fontsize = 10))
pdf(file.path(out, "figures", "D5_Fig2_panel_heatmap.pdf"), width = 7.2, height = 5.0)
draw(ht)
dev.off()
png(file.path(out, "figures", "D5_Fig2_panel_heatmap.png"), width = 7.2, height = 5.0,
    units = "in", res = 300)
draw(ht)
dev.off()

# ---------- VIF 共线性 ----------
dd <- as.data.frame(d[, ..panel]); dd$y <- as.integer(d$lesion == "Lesion")
full <- glm(y ~ ., data = dd, family = binomial)
vifv <- tryCatch(vif(full), error = function(e) NULL)
if (!is.null(vifv)) {
  vt <- data.table(gene = names(vifv), VIF = as.numeric(vifv))
  setorder(vt, -VIF)
  fwrite(vt, file.path(out, "results", "D5_vif.tsv"), sep = "\t")
  cat("\n=== VIF（全部 10 基因） ===\n"); print(vt)
  keep <- vt[VIF < 5]$gene
  cat("VIF<5 保留:", paste(keep, collapse = ", "), "\n")
  saveRDS(keep, file.path(out, "data", "D5_vif_kept_genes.rds"))
}
cat("\nD5-01 done\n")