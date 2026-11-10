# ============================================================
# D4-06 多因素 logistic 逐步回归交叉验证（计划书退路/best-practice 对照）
# 目的：对 TMSB4X 单基因模型与"三算法并集 top" 多基因模型做同一队列的
#       logistic 拟合 + AIC 逐步筛选 + ROC，确认 D4 收敛结果的稳健性
# 输出: results/D4_logistic_stepwise.tsv, D4_model_comparison.tsv
#       figures/D4_Fig7_logistic_roc.pdf/.png
# ============================================================
suppressPackageStartupMessages({
  library(Biobase); library(data.table); library(pROC); library(ggplot2)
})
out <- "D:/A生信分析/虚拟敲除/第4天执行结果"
es <- readRDS("D:/A生信分析/虚拟敲除/D1执行结果/D1_CNV_角膜新生血管_数据分析/data/processed/GSE241715_expression.rds")
x <- exprs(es); pd <- pData(es)
y <- as.integer(pd$group != "Control")

cand <- c("TMSB4X", "S100A13", "XIAP", "CAPG", "PRDX2", "VDAC1", "P4HB", "SOD1", "GJA1",
          "SERPINE1", "THBS1", "LOX")
cand <- intersect(cand, rownames(x))
d <- as.data.frame(t(x[cand, , drop = FALSE]))
d$y <- y
d <- cbind(d, t(x[c("TGFB1", "STAT3", "PXDN"), , drop = FALSE]) |> as.data.frame())

# 全变量多因素 logistic + AIC 双向逐步
full <- glm(y ~ ., data = d, family = binomial)
step_fit <- step(full, direction = "both", trace = 0)
co <- summary(step_fit)$coefficients
res <- data.table(term = rownames(co), estimate = co[, 1], se = co[, 2],
                  z = co[, 3], p = co[, 4], OR = exp(co[, 1]))
setorder(res, p)
fwrite(res, file.path(out, "results", "D4_logistic_stepwise.tsv"), sep = "\t")
cat("AIC stepwise logistic 保留变量:\n"); print(res)

# 模型比较：单基因 TMSB4X vs 逐步模型 vs 先验 SERPINE1
mk <- function(f, tag) {
  m <- glm(f, data = d, family = binomial)
  pr <- predict(m, type = "response")
  r <- roc(y, pr, quiet = TRUE, direction = "<")
  data.table(model = tag, n_terms = length(coef(m)) - 1,
             AIC = AIC(m),
             AUC = as.numeric(auc(r)),
             CI_low = as.numeric(ci.auc(r))[1], CI_high = as.numeric(ci.auc(r))[3],
             p_vs_null = anova(m, test = "Chisq")$`Pr(>Chi)`[2])
}
cmp <- rbind(
  mk(y ~ TMSB4X, "TMSB4X (single)"),
  mk(y ~ SERPINE1, "SERPINE1 (prior)"),
  mk(y ~ TMSB4X + SERPINE1, "TMSB4X + SERPINE1"),
  mk(as.formula(paste("y ~", paste(setdiff(res$term, "(Intercept)"), collapse = " + "))),
     "AIC stepwise (all candidates)"),
  mk(y ~ ., "full 15-gene model")
)
fwrite(cmp, file.path(out, "results", "D4_model_comparison.tsv"), sep = "\t")
print(cmp)

# ROC 曲线
rocs <- lapply(cmp$model, function(tag) {
  f <- switch(tag,
    "TMSB4X (single)" = y ~ TMSB4X,
    "SERPINE1 (prior)" = y ~ SERPINE1,
    "TMSB4X + SERPINE1" = y ~ TMSB4X + SERPINE1,
    "AIC stepwise (all candidates)" = as.formula(paste("y ~", paste(setdiff(res$term, "(Intercept)"), collapse = " + "))),
    y ~ .)
  m <- glm(f, data = d, family = binomial)
  r <- roc(y, predict(m, type = "response"), quiet = TRUE, direction = "<")
  data.table(model = tag, spec = rev(r$specificities), sens = rev(r$sensitivities))
})
roc_dt <- rbindlist(rocs)
p <- ggplot(roc_dt, aes(1 - spec, sens, colour = model)) +
  geom_abline(linetype = 3, colour = "grey70") +
  geom_line(linewidth = 0.9) +
  scale_colour_brewer(palette = "Set1") +
  coord_equal() +
  labs(title = "Training-set ROC of D4 candidate models (GSE241715, n = 21)",
       subtitle = "in-sample fit; D5 performs the prospective-style evaluation",
       x = "1 - Specificity", y = "Sensitivity", colour = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = c(0.72, 0.18))
ggsave(file.path(out, "figures", "D4_Fig7_logistic_roc.pdf"), p, width = 6.6, height = 5.6)
ggsave(file.path(out, "figures", "D4_Fig7_logistic_roc.png"), p, width = 6.6, height = 5.6, dpi = 300)
cat("D4-06 done\n")