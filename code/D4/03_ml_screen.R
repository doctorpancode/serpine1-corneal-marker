# ============================================================
# D4-03 三算法特征筛选：LASSO / RandomForest / SVM-RFE
# 输入: D3 candidate_pool.csv(167) + GSE241715 表达(log2 CPM) + 病变分组
# 输出: results/D4_lasso_coef.tsv / D4_rf_importance.tsv / D4_svmrfe_profile.tsv
#       results/D4_ml_selected.tsv / D4_ml_vote_matrix.tsv
#       data/D4_ml_selected.json
#       figures/D4_Fig3_lasso_cv / Fig4_rf_importance / Fig5_svmrfe / Fig6_ml_convergence
# ============================================================
suppressPackageStartupMessages({
  library(Biobase); library(data.table); library(glmnet)
  library(randomForest); library(caret); library(pROC); library(ggplot2)
})
set.seed(20261109)
out <- "D:/A生信分析/虚拟敲除/第4天执行结果"
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "results"), recursive = TRUE, showWarnings = FALSE)

es <- readRDS("D:/A生信分析/虚拟敲除/D1执行结果/D1_CNV_角膜新生血管_数据分析/data/processed/GSE241715_expression.rds")
x  <- exprs(es)
pd <- pData(es)
pool <- fread("D:/A生信分析/虚拟敲除/第3天执行结果/data/candidate_pool.csv")
genes <- intersect(pool$gene, rownames(x))
cat("候选池 167 -> 矩阵命中:", length(genes), "\n")

y <- factor(ifelse(pd$group == "Control", "Normal", "Lesion"), levels = c("Normal", "Lesion"))
tab <- table(y); print(tab)

X <- t(x[genes, , drop = FALSE])
# 去掉零方差特征
v <- apply(X, 2, var); keep <- names(v)[v > 1e-8]
if (length(keep) < ncol(X)) cat("零方差剔除:", ncol(X) - length(keep), "\n")
X <- X[, keep, drop = FALSE]
cs <- scale(X)                       # z-score（对二分类无信息泄漏风险，均值/方差全队列估计）
cs[is.na(cs)] <- 0

# ---------------- 1) LASSO ----------------
set.seed(20261109)
cvfit <- cv.glmnet(cs, y, family = "binomial", alpha = 1, nfolds = 10,
                   type.measure = "deviance", standardize = FALSE)
co <- as.matrix(coef(cvfit, s = "lambda.min"))
lasso <- data.table(gene = rownames(co), coef = as.numeric(co))
lasso <- lasso[gene != "(Intercept)"][order(-abs(coef))]
lasso[, selected := abs(coef) > 1e-8]
fwrite(lasso, file.path(out, "results", "D4_lasso_coef.tsv"), sep = "\t")
lasso_sel <- lasso[selected == TRUE]$gene
cat("\nLASSO lambda.min =", signif(cvfit$lambda.min, 5),
    " n_selected =", length(lasso_sel), ":", paste(lasso_sel, collapse = ", "), "\n")

# LASSO CV 图
pdf_dt <- data.table(lambda = cvfit$lambda, cvm = cvfit$cvm, cvsd = cvfit$cvsd,
                     nzero = cvfit$nzero)
p1 <- ggplot(pdf_dt, aes(log(lambda), cvm)) +
  geom_point(colour = "#4c72b0", size = 1.6) +
  geom_errorbar(aes(ymin = cvm - cvsd, ymax = cvm + cvsd), width = 0.02, colour = "grey55") +
  geom_vline(xintercept = log(cvfit$lambda.min), linetype = 2, colour = "#d62728") +
  geom_vline(xintercept = log(cvfit$lambda.1se), linetype = 3, colour = "grey30") +
  labs(title = "LASSO 10-fold CV (binomial deviance)",
       subtitle = paste0("lambda.min = ", signif(cvfit$lambda.min, 4),
                         "; nonzero = ", length(lasso_sel), " (dashed red)"),
       x = "log(lambda)", y = "Binomial deviance") + theme_bw(base_size = 11)
ggsave(file.path(out, "figures", "D4_Fig3_lasso_cv.pdf"), p1, width = 6.2, height = 4.6)
ggsave(file.path(out, "figures", "D4_Fig3_lasso_cv.png"), p1, width = 6.2, height = 4.6, dpi = 300)

# ---------------- 2) RandomForest ----------------
set.seed(20261109)
rf <- randomForest(x = cs, y = y, ntree = 500, importance = TRUE, proximity = FALSE)
imp <- data.table(gene = rownames(importance(rf)),
                  MDA = importance(rf)[, "MeanDecreaseAccuracy"],
                  MDG = importance(rf)[, "MeanDecreaseGini"])
imp <- imp[order(-MDG)]
fwrite(imp, file.path(out, "results", "D4_rf_importance.tsv"), sep = "\t")
rf_sel <- head(imp$gene, 20)
cat("RF top20 (MDG):", paste(rf_sel, collapse = ", "), "\n")
cat("RF OOB err =", signif(tail(rf$err.rate[, "OOB"], 1), 4), "\n")

p2 <- ggplot(head(imp, 20), aes(reorder(gene, MDG), MDG)) +
  geom_col(fill = "#55a868") + coord_flip() +
  labs(title = "RandomForest MeanDecreaseGini (top 20)", x = NULL, y = "MeanDecreaseGini") +
  theme_bw(base_size = 11)
ggsave(file.path(out, "figures", "D4_Fig4_rf_importance.pdf"), p2, width = 6.0, height = 5.4)
ggsave(file.path(out, "figures", "D4_Fig4_rf_importance.png"), p2, width = 6.0, height = 5.4, dpi = 300)

# ---------------- 3) SVM-RFE ----------------
ctrl <- rfeControl(functions = caretFuncs, method = "cv", number = 5,
                   repeats = 1, verbose = FALSE, allowParallel = FALSE)
subsets <- c(1, 2, 3, 5, 8, 12, 16, 20, 25, 30, length(rf_sel) + 5, min(60, ncol(cs)))
subsets <- sort(unique(subsets[subsets <= ncol(cs)]))
set.seed(20261109)
rfe_fit <- rfe(x = as.data.frame(cs), y = y, sizes = subsets,
               rfeControl = ctrl, method = "svmRadial",
               tuneLength = 3, metric = "ROC")
sv <- rfe_fit$variables
fwrite(as.data.table(sv), file.path(out, "results", "D4_svmrfe_profile.tsv"), sep = "\t")
best_n <- rfe_fit$bestSubset
svm_sel <- as.character(sv$var[sv$Variables == best_n])
cat("\nSVM-RFE best subset size =", best_n, " accuracy =", signif(max(rfe_fit$results$ROC), 4),
    "\n  genes:", paste(svm_sel, collapse = ", "), "\n")

prof <- as.data.table(rfe_fit$results)
p3 <- ggplot(prof, aes(Variables, ROC)) +
  geom_line(colour = "#4c72b0") + geom_point(colour = "#4c72b0", size = 2) +
  geom_point(data = prof[Variables == best_n], colour = "#d62728", size = 3.4) +
  labs(title = "SVM-RFE recursive feature elimination", subtitle = paste0("best subset = ", best_n, " genes"),
       x = "Number of features", y = "CV ROC (5-fold)") + theme_bw(base_size = 11)
ggsave(file.path(out, "figures", "D4_Fig5_svmrfe.pdf"), p3, width = 6.2, height = 4.6)
ggsave(file.path(out, "figures", "D4_Fig5_svmrfe.png"), p3, width = 6.2, height = 4.6, dpi = 300)

# ---------------- 4) 三算法投票收敛 + PPI hub 加权 ----------------
hub <- fread("D:/A生信分析/虚拟敲除/第4天执行结果/results/D4_hub_intersection.tsv")
vote <- data.table(gene = genes)
vote[, LASSO := as.integer(gene %in% lasso_sel)]
vote[, RF_top20 := as.integer(gene %in% rf_sel)]
vote[, SVM_RFE := as.integer(gene %in% svm_sel)]
vote[, PPI_hub := as.integer(gene %in% hub$gene)]
vote[, n_vote := LASSO + RF_top20 + SVM_RFE]
vote[, lasso_coef := lasso$coef[match(gene, lasso$gene)]]
vote[, rf_MDG := imp$MDG[match(gene, imp$gene)]]
vote <- vote[order(-n_vote, -PPI_hub, -rf_MDG)]
fwrite(vote, file.path(out, "results", "D4_ml_vote_matrix.tsv"), sep = "\t")
fwrite(vote[n_vote > 0], file.path(out, "results", "D4_ml_selected.tsv"), sep = "\t")

cat("\n=== 三算法投票结果（n_vote >= 2） ===\n")
print(vote[n_vote >= 2])

# 收敛：3 票 + PPI hub，或 2 票 + PPI hub，或 3 票全中
converged <- vote[n_vote == 3]
if (nrow(converged) == 0) converged <- vote[n_vote >= 2]
if (nrow(converged) == 0) converged <- vote[PPI_hub == 1]
# 最终 1-2 基因：按 n_vote、PPI_hub、|lasso_coef|、MDG 排序取前 2
final <- head(converged[order(-n_vote, -PPI_hub, -abs(lasso_coef), -rf_MDG)], 2)
if (nrow(final) == 0) final <- vote[order(-n_vote, -abs(lasso_coef))][1:2]
cat("\n>>> 最终收敛基因（1-2 个）:\n"); print(final)

# 强制记录先验三候选的落选/入选状态
prior <- vote[gene %in% c("SERPINE1", "THBS1", "LOX")]
cat("\n=== 先验候选基因在 D4 的位置 ===\n"); print(prior)

jsonlite::write_json(list(
  lasso_lambda_min = cvfit$lambda.min, lasso_selected = lasso_sel,
  rf_top20 = rf_sel, rf_oob_error = tail(rf$err.rate[, "OOB"], 1),
  svm_best_n = best_n, svm_selected = svm_sel,
  ppi_hub_intersection = hub$gene,
  final_genes = final$gene,
  prior_candidates = setNames(prior$n_vote, prior$gene)
), file.path(out, "data", "D4_ml_selected.json"), auto_unbox = TRUE, pretty = TRUE)

# ---------------- 5) 收敛图 ----------------
vm <- melt(vote[n_vote > 0][order(-n_vote)][1:min(30, .N)],
           id.vars = "gene", measure.vars = c("LASSO", "RF_top20", "SVM_RFE"))
vm[, gene := factor(gene, levels = rev(unique(vote[n_vote > 0][order(-n_vote)]$gene)))]
p4 <- ggplot(vm, aes(gene, variable, fill = factor(value))) +
  geom_tile(colour = "white") +
  scale_fill_manual(values = c("0" = "grey92", "1" = "#d62728"), guide = "none") +
  labs(title = "Three-algorithm evidence matrix", x = NULL, y = NULL,
       subtitle = paste0("final: ", paste(final$gene, collapse = " + "))) +
  theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave(file.path(out, "figures", "D4_Fig6_ml_convergence.pdf"), p4, width = 8.2, height = 4.2)
ggsave(file.path(out, "figures", "D4_Fig6_ml_convergence.png"), p4, width = 8.2, height = 4.2, dpi = 300)
cat("\nD4 ML done\n")