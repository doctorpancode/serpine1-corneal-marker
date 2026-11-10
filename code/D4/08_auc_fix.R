# ============================================================
# D4-08 修正 AUC 估计（分离型数据的 DeLong CI）并重绘 ROC
# 背景：TMSB4X 在 GSE241715 中 7 例 Control 全部高于 14 例 Lesion（完全分离）
#       → logistic 系数发散，AUC=1.0，DeLong 方差为 0（CI 退化为 1-1）
#       论文口径必须写成 "90% CI (bootstrap)" + "training-set, upper bound"
# ============================================================
suppressPackageStartupMessages({
  library(Biobase); library(data.table); library(pROC); library(ggplot2)
})
set.seed(20261109)
out <- "D:/A生信分析/虚拟敲除/第4天执行结果"
es <- readRDS("D:/A生信分析/虚拟敲除/D1执行结果/D1_CNV_角膜新生血管_数据分析/data/processed/GSE241715_expression.rds")
x <- exprs(es); pd <- pData(es)
y <- as.integer(pd$group != "Control")
genes <- c("TMSB4X", "SERPINE1", "S100A13", "XIAP", "CAPG", "PRDX2", "VDAC1", "P4HB", "SOD1", "GJA1")
d <- data.table(t(x[genes, , drop = FALSE])); d[, y := y]

models <- list(
  "TMSB4X (single)"  = "y ~ TMSB4X",
  "SERPINE1 (prior)" = "y ~ SERPINE1",
  "S100A13"          = "y ~ S100A13",
  "TMSB4X + SERPINE1" = "y ~ TMSB4X + SERPINE1",
  "10-gene panel"    = paste("y ~", paste(genes, collapse = " + "))
)

boot_auc <- function(score, yy, B = 2000) {
  n <- length(yy)
  v <- numeric(0)
  for (i in 1:B) {
    idx <- sample(n, n, TRUE)
    if (length(unique(yy[idx])) < 2) next
    r <- pROC::roc(yy[idx], score[idx], quiet = TRUE, direction = "<")
    v <- c(v, as.numeric(pROC::auc(r)))
  }
  quantile(v, c(0.025, 0.975), na.rm = TRUE)
}

rows <- list(); roc_dt <- list()
for (tag in names(models)) {
  f <- as.formula(models[[tag]])
  rank_trans <- FALSE
  if (tag %in% c("TMSB4X (single)", "SERPINE1 (prior)", "S100A13")) {
    g <- all.vars(f)[2]
    score <- d[[g]] * if (g == "TMSB4X") -1 else 1   # 方向按生物学预期固定
    m <- NULL
  } else {
    m <- glm(f, data = d, family = binomial)
    score <- fitted(m)
  }
  r <- pROC::roc(y, score, quiet = TRUE, direction = "<")
  ci_delong <- as.numeric(pROC::ci.auc(r, method = "delong"))
  ci_boot <- boot_auc(score, y)
  rows[[tag]] <- data.table(
    model = tag,
    n_terms = if (is.null(m)) 1 else length(coef(m)) - 1,
    AIC = if (is.null(m)) NA_real_ else AIC(m),
    AUC = as.numeric(pROC::auc(r)),
    CI_delong_low = ci_delong[1], CI_delong_high = ci_delong[3],
    CI_boot_low = as.numeric(ci_boot[1]), CI_boot_high = as.numeric(ci_boot[2])
  )
  roc_dt[[tag]] <- data.table(model = tag, spec = rev(r$specificities), sens = rev(r$sensitivities))
}
cmp <- rbindlist(rows)
fwrite(cmp, file.path(out, "results", "D4_model_comparison.tsv"), sep = "\t")
print(cmp)

p <- ggplot(rbindlist(roc_dt), aes(1 - spec, sens, colour = model)) +
  geom_abline(linetype = 3, colour = "grey70") +
  geom_line(linewidth = 0.95) + coord_equal() +
  scale_colour_brewer(palette = "Set1") +
  labs(title = "Training-set ROC (GSE241715, n = 21: 7 controls / 14 lesions)",
       subtitle = "TMSB4X fully separates groups -> AUC 1.00 is an upper bound; bootstrap 95% CI reported in table",
       x = "1 - Specificity", y = "Sensitivity", colour = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "inside", legend.position.inside = c(0.72, 0.2))
ggsave(file.path(out, "figures", "D4_Fig7_logistic_roc.pdf"), p, width = 6.8, height = 5.8)
ggsave(file.path(out, "figures", "D4_Fig7_logistic_roc.png"), p, width = 6.8, height = 5.8, dpi = 300)
cat("D4-08 done\n")