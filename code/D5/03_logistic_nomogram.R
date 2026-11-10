# ============================================================
# D5-03 诊断模型：多因素 logistic + ROC + 校准曲线 + 列线图 + 决策曲线(DCA)
# 训练集：GSE241715 (n=21, 7 Control / 14 Lesion)
# 外部验证：GSE204791 基质室 (n=15, 7 Control / 8 KC)
# 注意：n=21 时任何 3 基因以上的模型都会过拟合 → 以 1-2 基因模型为主，报告训练集上界
# 输出：results/D5_model_*.tsv, data/D5_models.rds
#       figures/D5_Fig4_ROC.pdf/.png, D5_Fig5_calibration.pdf/.png, D5_Fig6_nomogram.pdf/.png
# ============================================================
suppressPackageStartupMessages({
  library(Biobase); library(data.table); library(ggplot2); library(pROC)
  library(rms); library(patchwork)
})
set.seed(20261109)
out <- "D:/A生信分析/虚拟敲除/第5天执行结果"
DATA <- "D:/A生信分析/虚拟敲除/D1执行结果/D1_CNV_角膜新生血管_数据分析/data/processed"

# ---------------- 训练集 ----------------
es <- readRDS(file.path(DATA, "GSE241715_expression.rds"))
x <- exprs(es); pd <- pData(es)
panel <- c("TMSB4X", "SERPINE1", "S100A13", "XIAP", "CAPG", "PRDX2", "VDAC1", "P4HB", "SOD1", "GJA1")
panel <- intersect(panel, rownames(x))
tr <- data.table(t(x[panel, , drop = FALSE]))
tr[, y := as.integer(pd$group != "Control")]
tr[, group := pd$group]

models <- list(
  M1_TMSB4X     = c("TMSB4X"),
  M2_SERPINE1   = c("SERPINE1"),
  M3_TwoCore    = c("TMSB4X", "SERPINE1"),
  M4_ThreeTop   = c("TMSB4X", "SERPINE1", "S100A13"),
  M5_Panel10    = panel
)

rocl <- list(); fitstats <- list()
rocobjs <- list()
for (nm in names(models)) {
  g <- models[[nm]]
  f <- as.formula(paste("y ~", paste(g, collapse = " + ")))
  m <- glm(f, data = tr, family = binomial)
  pr <- predict(m, type = "response")
  r <- roc(tr$y, pr, quiet = TRUE)
  ci <- tryCatch(as.numeric(ci.auc(r, method = "delong")), error = function(e) c(NA, NA, NA))
  ll <- as.numeric(logLik(m))
  k <- length(coef(m))
  aic <- -2 * ll + 2 * k
  bic <- -2 * ll + k * log(nrow(tr))
  # Hosmer-Lemeshow
  hl <- tryCatch({
    q <- unique(quantile(pr, probs = seq(0, 1, 0.2), na.rm = TRUE))
    gq <- cut(pr, breaks = q, include.lowest = TRUE, labels = FALSE)
    obs1 <- tapply(tr$y, gq, sum); exp1 <- tapply(pr, gq, sum); n1 <- tapply(tr$y, gq, length)
    chis <- sum((obs1 - exp1)^2 / pmax(exp1, 1e-9))
    data.table(HL_chisq = chis, HL_df = length(obs1) - 2,
               HL_p = pchisq(chis, length(obs1) - 2, lower.tail = FALSE))
  }, error = function(e) data.table(HL_chisq = NA, HL_df = NA, HL_p = NA))
  rocl[[nm]] <- data.table(model = nm, n_gene = length(g), AUC = as.numeric(auc(r)),
                           CI_low = ci[1], CI_high = ci[3], AIC = aic, BIC = bic, HL_p = hl$HL_p)
  fitstats[[nm]] <- data.table(model = nm, term = names(coef(m)),
                               beta = as.numeric(coef(m)),
                               OR = exp(as.numeric(coef(m))),
                               p = summary(m)$coefficients[, 4])
  rocobjs[[nm]] <- r
}
R <- rbindlist(rocl); F <- rbindlist(fitstats)
fwrite(R, file.path(out, "results", "D5_model_ROCAUC.tsv"), sep = "\t")
fwrite(F, file.path(out, "results", "D5_model_coefficients.tsv"), sep = "\t")
print(R[order(-AUC)]); cat("\n"); print(F)

# ---------------- 主模型（双主干）系数 ----------------
main <- glm(y ~ TMSB4X + SERPINE1, data = tr, family = binomial)
cat("\n=== 主模型 y ~ TMSB4X + SERPINE1 ===\n"); print(summary(main)$coefficients)
print(cbind(OR = exp(coef(main)), exp(confint(main))))
fwrite(data.table(term = names(coef(main)), beta = coef(main), OR = exp(coef(main)),
                  OR_low = exp(confint(main))[, 1], OR_high = exp(confint(main))[, 2],
                  p = summary(main)$coefficients[, 4]),
       file.path(out, "results", "D5_main_model_OR.tsv"), sep = "\t")

# 单因素（每个基因单独 logistic：OR 可解释、无共线性）
uni <- rbindlist(lapply(panel, function(g) {
  m <- glm(tr[[g]] * 0 + tr$y ~ tr[[g]], family = binomial)
  s <- summary(m)$coefficients
  data.table(gene = g, beta = s[2, 1], OR = exp(s[2, 1]),
             OR_low = exp(confint(m))[2, 1], OR_high = exp(confint(m))[2, 2],
             p = s[2, 4])
}))
setorder(uni, p)
fwrite(uni, file.path(out, "results", "D5_univariable_OR.tsv"), sep = "\t")
print(uni)

# ---------------- 决策曲线分析（手工实现） ----------------
dca <- rbindlist(lapply(names(rocobjs), function(nm) {
  g <- models[[nm]]
  pr <- predict(glm(as.formula(paste("y ~", paste(g, collapse = " + "))), data = tr,
                    family = binomial), type = "response")
  th <- seq(0.01, 0.99, 0.01)
  tms <- data.table(model = nm, threshold = th,
                    nb_model = sapply(th, function(t) sum(pr >= t & tr$y == 1) / nrow(tr) -
                                        sum(pr >= t & tr$y == 0) / nrow(tr) * t / (1 - t)))
  tms
}))
dca <- rbind(dca, data.table(model = "Treat_all", threshold = seq(0.01, 0.99, 0.01),
                             nb_model = sapply(seq(0.01, 0.99, 0.01), function(t)
                               mean(tr$y) - (1 - mean(tr$y)) * t / (1 - t))))
fwrite(dca, file.path(out, "results", "D5_DCA.tsv"), sep = "\t")

# ---------------- 图4：ROC ----------------
rocdf <- rbindlist(lapply(names(rocobjs), function(nm) {
  r <- rocobjs[[nm]]
  data.table(model = nm, fpr = 1 - r$specificities, tpr = r$sensitivities)
}))
lab <- R[, .(model, lab = sprintf("%s (AUC %.3f)", model, AUC))]
rocdf <- merge(rocdf, lab, by = "model")
pR <- ggplot(rocdf, aes(fpr, tpr, colour = lab)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
  geom_line(linewidth = 0.85, alpha = 0.9) +
  scale_linetype_manual(values = c(1,2,3,4,6)) +
  scale_colour_brewer(palette = "Set1") +
  labs(title = "Diagnostic ROC, GSE241715 training set (n = 21)",
       subtitle = "AUC values are training-set upper bounds; no held-out split at n = 21",
       x = "1 - Specificity", y = "Sensitivity", colour = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = c(0.68, 0.25),
                                   legend.background = element_rect(fill = alpha("white", .8)))
ggsave(file.path(out, "figures", "D5_Fig4_ROC.pdf"), pR, width = 6.6, height = 5.4)
ggsave(file.path(out, "figures", "D5_Fig4_ROC.png"), pR, width = 6.6, height = 5.4, dpi = 300)

# ---------------- 图5：校准曲线（bootstrap 平滑） ----------------
calib <- rbindlist(lapply(c("M1_TMSB4X", "M2_SERPINE1", "M3_TwoCore"), function(nm) {
  g <- models[[nm]]
  pr <- predict(glm(as.formula(paste("y ~", paste(g, collapse = " + "))), data = tr,
                    family = binomial), type = "response")
  data.table(model = nm, pred = pr, obs = tr$y)
}))
pC <- ggplot(calib, aes(pred, obs, colour = model)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_smooth(method = "loess", se = FALSE, span = 1.1, linewidth = 0.9) +
  geom_point(size = 1.4, alpha = 0.75) +
  coord_cartesian(xlim = c(0, 1), ylim = c(-0.05, 1.05)) +
  scale_colour_brewer(palette = "Dark2") +
  labs(title = "Calibration of diagnostic models (GSE241715)",
       subtitle = "Observed vs predicted lesion probability; loess smoother",
       x = "Predicted probability", y = "Observed lesion status") +
  theme_bw(base_size = 10)
ggsave(file.path(out, "figures", "D5_Fig5_calibration.pdf"), pC, width = 6.4, height = 5.2)
ggsave(file.path(out, "figures", "D5_Fig5_calibration.png"), pC, width = 6.4, height = 5.2, dpi = 300)

# ---------------- 图6：列线图 ----------------
dd <- as.data.frame(tr)
dd$y <- factor(dd$y, levels = c(0, 1), labels = c("Normal", "Lesion"))
ddist <- datadist(dd); options(datadist = "ddist")
lrm1 <- lrm(y ~ TMSB4X + SERPINE1, data = dd, x = TRUE, y = TRUE)
cat("\nlrm C-index:", lrm1$stats["C"], "\n")
pdf(file.path(out, "figures", "D5_Fig6_nomogram.pdf"), width = 7.6, height = 3.9)
plot(nomogram(lrm1, fun = plogis, funlabel = "Lesion probability",
              fun.at = c(0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99)),
     xfrac = 0.32, cex.axis = 0.8, cex.var = 0.85)
dev.off()
png(file.path(out, "figures", "D5_Fig6_nomogram.png"), width = 7.6, height = 3.9, units = "in", res = 300)
plot(nomogram(lrm1, fun = plogis, funlabel = "Lesion probability",
              fun.at = c(0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99)),
     xfrac = 0.32, cex.axis = 0.8, cex.var = 0.85)
dev.off()
# 列线图评分表（便于引用具体分值）
png(file.path(out, "figures", "D5_Fig6B_nomogram_influence.pdf"), width = 7.6, height = 3.4,
    units = "in", res = 300)
plot(nomogram(lrm1, fun = plogis), xfrac = 0.32, cex.axis = 0.75, cex.var = 0.8)
dev.off()

# ---------------- 图7：DCA ----------------
pD <- ggplot(dca, aes(threshold, nb_model, colour = model)) +
  geom_line(linewidth = 0.8) +
  coord_cartesian(xlim = c(0, 0.9), ylim = c(-0.1, 0.65)) +
  scale_colour_brewer(palette = "Set1") +
  labs(title = "Decision curve analysis (GSE241715)",
       subtitle = "Net benefit of the diagnostic models across threshold probabilities",
       x = "Threshold probability", y = "Net benefit", colour = NULL) +
  theme_bw(base_size = 10)
ggsave(file.path(out, "figures", "D5_Fig7_DCA.pdf"), pD, width = 6.6, height = 5.0)
ggsave(file.path(out, "figures", "D5_Fig7_DCA.png"), pD, width = 6.6, height = 5.0, dpi = 300)

saveRDS(list(models = models, main = main, lrm = lrm1, roc = rocobjs),
        file.path(out, "data", "D5_models.rds"))
cat("\nD5-03 done\n")