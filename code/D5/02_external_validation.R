# ============================================================
# D5-02 独立检验队列：GSE204791 上皮/基质（人眼角膜）
# 背景：GSE241715 仅 21 例，模型 AUC 达 1.0 属训练集上界，必须在独立队列验证
# 输出：results/D5_GSE204791_{epi,stroma}_stats.tsv
#       results/D5_GSE204791_ROC.tsv
#       figures/D5_Fig3_external_validation.pdf/.png
# ============================================================
suppressPackageStartupMessages({
  library(Biobase); library(data.table); library(ggplot2); library(pROC)
})
out <- "D:/A生信分析/虚拟敲除/第5天执行结果"
DATA <- "D:/A生信分析/虚拟敲除/D1执行结果/D1_CNV_角膜新生血管_数据分析/data/processed"

es <- readRDS(file.path(DATA, "GSE204791_expression.rds"))
ann <- read.delim("D:/A生信分析/虚拟敲除/D1执行结果/D1_CNV_角膜新生血管_数据分析/results/GSE204791_probe_annotation.tsv",
                  stringsAsFactors = FALSE)
ss <- fread(file.path(DATA, "GSE204791_sample_sheet.tsv"))
x <- exprs(es)
# 探针 -> 基因：取每个基因最大平均表达探针
ann <- ann[ann$SYMBOL != "" & !is.na(ann$SYMBOL), ]
ann <- ann[ann$probe %in% rownames(x), ]
mx <- rowMeans(x[ann$probe, , drop = FALSE])
ann$mu <- mx
ann <- ann[order(ann$SYMBOL, -ann$mu), ]
ann <- ann[!duplicated(ann$SYMBOL), ]
xg <- x[ann$probe, , drop = FALSE]
rownames(xg) <- ann$SYMBOL
cat("基因矩阵:", dim(xg), "\n")

panel <- c("TMSB4X", "SERPINE1", "S100A13", "XIAP", "CAPG", "PRDX2", "VDAC1", "P4HB", "SOD1", "GJA1")
panel <- intersect(panel, rownames(xg))
miss <- setdiff(c("TMSB4X", "SERPINE1"), rownames(xg))
if (length(miss)) cat("缺失:", miss, "\n")

ss <- ss[match(colnames(xg), ss$gsm), ]
stopifnot(!any(is.na(ss$gsm)))
d <- data.table(t(xg[panel, , drop = FALSE]))
d[, gsm := colnames(xg)]
d[, compartment := ss$compartment]
d[, group := ss$group]

res <- list(); rocres <- list()
for (cmp in unique(ss$compartment)) {
  dd <- d[compartment == cmp]
  dd[, y := as.integer(group == "KC")]
  cat("\n==== ", cmp, ": n=", nrow(dd), " (KC=", sum(dd$y), ", Ctrl=", sum(dd$y == 0), ") ====\n", sep = "")
  rr <- rbindlist(lapply(panel, function(g) {
    v <- dd[[g]]
    p <- suppressWarnings(wilcox.test(v ~ dd$y)$p.value)
    data.table(compartment = cmp, gene = g,
               mean_ctrl = mean(v[dd$y == 0]), mean_kc = mean(v[dd$y == 1]),
               delta = mean(v[dd$y == 1]) - mean(v[dd$y == 0]),
               p_wilcox = p)
  }))
  rr[, adjP := p.adjust(p_wilcox, "BH")]
  setorder(rr, p_wilcox)
  res[[cmp]] <- rr
  print(rr[, .(gene, delta = round(delta, 3), p = signif(p_wilcox, 3), adjP = signif(adjP, 3))])

  for (g in panel) {
    r <- tryCatch(roc(dd$y, dd[[g]], quiet = TRUE, direction = "<"), error = function(e) NULL)
    if (is.null(r)) next
    if (auc(r) < 0.5) r <- roc(dd$y, dd[[g]], quiet = TRUE, direction = ">")
    ci <- tryCatch(as.numeric(ci.auc(r, method = "delong")), error = function(e) c(NA, NA, NA))
    rocres[[paste(cmp, g)]] <- data.table(compartment = cmp, gene = g,
                                          AUC = as.numeric(auc(r)),
                                          CI_low = ci[1], CI_high = ci[3])
  }
}
RS <- rbindlist(res); rocr <- rbindlist(rocres)
fwrite(RS, file.path(out, "results", "D5_GSE204791_stats.tsv"), sep = "\t")
fwrite(rocr, file.path(out, "results", "D5_GSE204791_ROC.tsv"), sep = "\t")
print(rocr[order(-AUC)])

# ---- 图：外部验证（小提琴 + AUC 条形） ----
pl <- melt(d[group %in% c("Control", "KC")],
           id.vars = c("gsm", "compartment", "group"), measure.vars = panel,
           variable.name = "gene", value.name = "expr")
pl[, grp := factor(ifelse(group == "KC", "Keratoconus", "Control"),
                   levels = c("Control", "Keratoconus"))]
pA <- ggplot(pl, aes(grp, expr, fill = grp)) +
  geom_violin(scale = "width", alpha = 0.55, colour = NA) +
  geom_boxplot(width = 0.13, outlier.shape = NA, fill = "white", alpha = 0.9) +
  geom_jitter(width = 0.08, size = 1, alpha = 0.65) +
  facet_grid(gene ~ compartment, scales = "free_y") +
  scale_fill_manual(values = c(Control = "#4c72b0", Keratoconus = "#d62728")) +
  labs(title = "External validation set GSE204791 (human cornea, KC vs Control)",
       subtitle = "per-compartment Wilcoxon; n = 8 vs 8 epithelium, 7 vs 8 stroma",
       x = NULL, y = "log2 signal") +
  theme_bw(base_size = 9) + theme(legend.position = "none",
                                  axis.text.x = element_text(angle = 30, hjust = 1))
pB <- ggplot(rocr, aes(reorder(gene, AUC), AUC, fill = compartment)) +
  geom_col(position = position_dodge(0.7), width = 0.65) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), position = position_dodge(0.7), width = 0.2) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey40") +
  coord_flip() + ylim(0, 1) +
  scale_fill_manual(values = c(epithelium = "#4c72b0", stroma = "#55a868")) +
  labs(title = "Discriminative ability in GSE204791", x = NULL, y = "AUC (95% DeLong CI)") +
  theme_bw(base_size = 9)
library(patchwork)
p <- wrap_plots(pA, pB, ncol = 2, widths = c(1.35, 1))
ggsave(file.path(out, "figures", "D5_Fig3_external_validation.pdf"), p, width = 12.5, height = 6.6)
ggsave(file.path(out, "figures", "D5_Fig3_external_validation.png"), p, width = 12.5, height = 6.6, dpi = 300)
cat("\nD5-02 done\n")