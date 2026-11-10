# ============================================================
# D4-09b 跨队列核验（修正版）
# 修 1: long3 的熔解变量名与列名撞车 -> 用 measure() 显式指定
# 修 2: GSE204791 是 quantile-normalised log2 signal（2^signal 才是线性强度），
#       跨分层比较必须在 log2 空间做，不能直接报原始 signal 均值
# ============================================================
suppressPackageStartupMessages({
  library(Biobase); library(data.table); library(ggplot2); library(patchwork)
})
out <- "D:/A生信分析/虚拟敲除/第4天执行结果"
genes <- c("TMSB4X", "SERPINE1", "S100A13", "XIAP", "CAPG", "PRDX2", "VDAC1", "P4HB", "SOD1", "GJA1")

# ---------- A) GSE204791 ----------
e2 <- readRDS("D:/A生信分析/虚拟敲除/D1执行结果/D1_CNV_角膜新生血管_数据分析/data/processed/GSE204791_expression.rds")
ann <- fread("D:/A生信分析/虚拟敲除/D1执行结果/D1_CNV_角膜新生血管_数据分析/results/GSE204791_probe_annotation.tsv")
ss  <- fread("D:/A生信分析/虚拟敲除/D1执行结果/D1_CNV_角膜新生血管_数据分析/data/processed/GSE204791_sample_sheet.tsv")
ex <- exprs(e2)
m <- merge(data.table(probe = rownames(ex), mean_sig = rowMeans(ex)), ann, by = "probe")
m <- m[SYMBOL %in% genes][order(SYMBOL, -mean_sig)]
m <- m[!duplicated(SYMBOL)]
sub <- ex[m$probe, , drop = FALSE]
rownames(sub) <- m$SYMBOL
d2 <- as.data.table(t(sub)); d2[, gsm := colnames(sub)]
d2 <- merge(d2, ss[, .(gsm, compartment, group)], by = "gsm")
d2[, compartment := factor(compartment, levels = c("Epithelium", "Stroma"))]
long2 <- melt(d2, id.vars = c("gsm", "compartment", "group"),
              measure.vars = genes, variable.name = "gene", value.name = "expr")
r2 <- long2[, {
  a <- expr[group == "KC"]; b <- expr[group == "Control"]
  .(n_KC = length(a), n_Ctrl = length(b), mean_KC = mean(a), mean_Ctrl = mean(b),
    delta = mean(a) - mean(b),
    p = suppressWarnings(tryCatch(wilcox.test(a, b)$p.value, error = function(e) NA_real_)))
}, by = .(compartment, gene)]
r2[, cohort := "GSE204791 (human cornea, log2 signal)"]
setorder(r2, compartment, gene)
print(r2)

g2 <- ggplot(long2, aes(group, expr, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
  geom_jitter(width = 0.13, size = 0.75, alpha = 0.6) +
  facet_wrap(~ compartment + gene, ncol = 5, scales = "free_y") +
  scale_fill_manual(values = c(Control = "#4c72b0", KC = "#d62728")) +
  labs(title = "A  GSE204791 (human cornea): D4 panel, KC vs Control",
       subtitle = "epithelium 8/8, stroma 8/7 (KC/Control); quantile-normalised log2 signal",
       x = NULL, y = "log2 signal") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none", strip.text = element_text(size = 6.5),
        axis.text.x = element_text(angle = 45, hjust = 1))

# ---------- B) GSE89537 ----------
e3 <- readRDS("D:/A生信分析/虚拟敲除/D1执行结果/D1_CNV_角膜新生血管_数据分析/data/processed/GSE89537_expression.rds")
ex3 <- exprs(e3); pd3 <- pData(e3)
g3n <- intersect(genes, rownames(ex3))
d3 <- data.table(t(ex3[g3n, , drop = FALSE]))
d3[, timepoint := factor(pd3$group, levels = c("Control_day0", "CNV_day5", "CNV_day15"))]
long3 <- melt(d3, id.vars = "timepoint", measure.vars = g3n,
              variable.name = "gene", value.name = "expr")
r3 <- long3[, .(n = .N, mean_expr = mean(expr)), by = .(gene, timepoint)]
w3 <- dcast(r3, gene ~ timepoint, value.var = "mean_expr")
w3[, log2FC_d5 := `CNV_day5` - Control_day0]
w3[, log2FC_d15 := `CNV_day15` - Control_day0]
w3[, cohort := "GSE89537 (mouse alkali-burn CNV, pooled)"]
print(w3)

g3 <- ggplot(long3, aes(timepoint, expr, group = gene, colour = gene)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.4) +
  labs(title = "B  GSE89537 (mouse): same panel across the CNV time course",
       subtitle = "pooled group means (2 biological replicates per point)",
       x = NULL, y = "log2 expression") +
  theme_bw(base_size = 9) +
  theme(legend.position = "right", legend.text = element_text(size = 7),
        legend.title = element_blank(), axis.text.x = element_text(angle = 20, hjust = 1))

pA <- g2 / g3 + plot_layout(heights = c(1.15, 1))
fig <- wrap_plots(g2, g3, ncol = 1, heights = c(1.15, 1))
ggsave(file.path(out, "figures", "D4_Fig8_crosscohort_final.pdf"), fig, width = 12, height = 9.5)
ggsave(file.path(out, "figures", "D4_Fig8_crosscohort_final.png"), fig, width = 12, height = 9.5, dpi = 300)

fwrite(r2, file.path(out, "results", "D4_crosscohort_final.tsv"), sep = "\t")
fwrite(w3, file.path(out, "results", "D4_crosscohort_mouse_final.tsv"), sep = "\t")
cat("\nD4-09b done\n")