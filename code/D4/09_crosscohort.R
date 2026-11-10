# ============================================================
# D4-09 跨队列方向一致性核验（D4 收尾的诚实性检查）
# 为什么需要：GSE241715 中 TMSB4X 完全分离（AUC 1.0），但那是 21 例小队列的
#             阈值效应。若换队列方向反转，就不能作为 CNV 的主候选。
# 检查对象: GSE204791（人角膜上皮/基质 KC vs Control）、GSE89537（鼠碱烧伤 CNV 全角膜 day0/day5/day15）
# 输出: results/D4_crosscohort_final_gene.tsv + figures/D4_Fig8_crosscohort_final.pdf/.png
# ============================================================
suppressPackageStartupMessages({
  library(Biobase); library(data.table); library(ggplot2); library(patchwork)
})
out <- "D:/A生信分析/虚拟敲除/第4天执行结果"
genes <- c("TMSB4X", "SERPINE1", "S100A13", "XIAP", "CAPG", "PRDX2", "VDAC1", "P4HB", "SOD1", "GJA1")

# ---------- A) GSE204791 (人, 上皮/基质) ----------
e2 <- readRDS("D:/A生信分析/虚拟敲除/D1执行结果/D1_CNV_角膜新生血管_数据分析/data/processed/GSE204791_expression.rds")
ann <- fread("D:/A生信分析/虚拟敲除/D1执行结果/D1_CNV_角膜新生血管_数据分析/results/GSE204791_probe_annotation.tsv")
ss <- fread("D:/A生信分析/虚拟敲除/D1执行结果/D1_CNV_角膜新生血管_数据分析/data/processed/GSE204791_sample_sheet.tsv")
ex <- exprs(e2)
# 探针 -> 基因（同基因多探针取最大平均信号）
m <- merge(data.table(probe = rownames(ex), mean_iso = rowMeans(ex)), ann, by = "probe")
m <- m[SYMBOL %in% genes][order(SYMBOL, -mean_iso)]
m <- m[!duplicated(SYMBOL)]
sub <- ex[m$probe, , drop = FALSE]
rownames(sub) <- m$SYMBOL
d2 <- as.data.table(t(sub))
d2[, gsm := rownames(t(sub))]
d2 <- merge(d2, ss[, .(gsm, compartment, group)], by = "gsm")
d2[, compartment := factor(compartment, levels = c("Epithelium", "Stroma"))]
long2 <- melt(d2, id.vars = c("gsm", "compartment", "group"), variable.name = "gene", value.name = "expr")
r2 <- long2[, {
  a <- expr[group == "KC"]; b <- expr[group == "Control"]
  .(n_KC = length(a), n_Ctrl = length(b),
    mean_KC = mean(a), mean_Ctrl = mean(b), delta = mean(a) - mean(b),
    p = tryCatch(wilcox.test(a, b)$p.value, error = function(e) NA_real_))
}, by = .(compartment, gene)]
r2[, cohort := "GSE204791 (human cornea)"]
print(r2)

g2 <- ggplot(long2, aes(group, expr, fill = group)) +
  geom_boxplot(outlier.size = 0.6, width = 0.6, alpha = 0.85) +
  geom_jitter(width = 0.12, size = 0.7, alpha = 0.55) +
  facet_grid(gene ~ compartment, scales = "free_y") +
  scale_fill_manual(values = c(Control = "#4c72b0", KC = "#d62728")) +
  labs(title = "A  GSE204791 (human): final-gene panel, KC vs Control",
       subtitle = "epithelium n=8/8, stroma n=8/7 (KC/Control); log2 quantile-normalised signal",
       x = NULL, y = "log2 expression") +
  theme_bw(base_size = 9) + theme(legend.position = "none",
                                  strip.text.y = element_text(size = 7, angle = 0))

# ---------- B) GSE89537 (鼠, 碱烧伤 CNV 时序) ----------
e3 <- readRDS("D:/A生信分析/虚拟敲除/D1执行结果/D1_CNV_角膜新生血管_数据分析/data/processed/GSE89537_expression.rds")
ex3 <- exprs(e3); pd3 <- pData(e3)
g3n <- intersect(genes, rownames(ex3))
d3 <- data.table(t(ex3[g3n, , drop = FALSE]))
d3[, timepoint := factor(pd3$group, levels = c("Control_day0", "CNV_day5", "CNV_day15"))]
long3 <- melt(d3, id.vars = "timepoint", variable.name = "gene", value.name = "expr")
r3 <- long3[, .(n = .N, mean_expr = mean(expr), sd = sd(expr)), by = .(gene, timepoint)]
w3 <- dcast(r3, gene ~ timepoint, value.var = "mean_expr")
setDT(w3)
w3[, log2FC_d5_vs_d0 := `CNV_day5` - Control_day0]
w3[, log2FC_d15_vs_d0 := `CNV_day15` - Control_day0]
w3[, cohort := "GSE89537 (mouse alkali-burn CNV, pooled n=2/group)"]
print(w3)

g3 <- ggplot(long3, aes(timepoint, expr, group = gene, colour = gene)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.2) +
  labs(title = "B  GSE89537 (mouse): same panel over the CNV time course",
       subtitle = "pooled group means (2 biological replicates per point; no replicates -> descriptive only)",
       x = NULL, y = "log2 expression") +
  theme_bw(base_size = 9) + theme(legend.position = "right",
                                  legend.text = element_text(size = 7),
                                  legend.title = element_blank())

fig <- wrap_plots(g2, g3, ncol = 2, widths = c(1.35, 1))
ggsave(file.path(out, "figures", "D4_Fig8_crosscohort_final.pdf"), fig, width = 13.5, height = 7.2)
ggsave(file.path(out, "figures", "D4_Fig8_crosscohort_final.png"), fig, width = 13.5, height = 7.2, dpi = 300)

fwrite(r2, file.path(out, "results", "D4_crosscohort_final.tsv"), sep = "\t")
fwrite(w3, file.path(out, "results", "D4_crosscohort_mouse_final.tsv"), sep = "\t")
cat("\nD4-09 done\n")