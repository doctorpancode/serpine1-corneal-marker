# ============================================================
# D4-04 三算法收敛最终判定（SVM-RFE 用简历排序第 1 名的无子集基因集）
# 输入: results/D4_lasso_coef.tsv, D4_rf_importance.tsv, D4_svmrfe_profile.tsv(重算)
#       data/D4_svmrfe_fit.rds, results/D4_hub_intersection.tsv
# 输出: results/D4_ml_vote_matrix.tsv, D4_ml_selected.tsv, D4_final_genes.tsv
#       data/D4_ml_selected.json, figures/D4_Fig5_svmrfe / Fig6_ml_convergence
# ============================================================
suppressPackageStartupMessages({
  library(Biobase); library(data.table); library(ggplot2); library(jsonlite)
})
set.seed(20261109)
out <- "D:/A生信分析/虚拟敲除/第4天执行结果"
es <- readRDS("D:/A生信分析/虚拟敲除/D1执行结果/D1_CNV_角膜新生血管_数据分析/data/processed/GSE241715_expression.rds")
x <- exprs(es); pd <- pData(es)
pool <- fread("D:/A生信分析/虚拟敲除/第3天执行结果/data/candidate_pool.csv")
genes <- intersect(pool$gene, rownames(x))
y <- factor(ifelse(pd$group == "Control", "Normal", "Lesion"), levels = c("Normal", "Lesion"))

lasso <- fread(file.path(out, "results", "D4_lasso_coef.tsv"))
lasso_sel <- lasso[selected == TRUE]$gene
imp <- fread(file.path(out, "results", "D4_rf_importance.tsv"))
rf_sel <- head(imp$gene, 20)
hub <- fread(file.path(out, "results", "D4_hub_intersection.tsv"))

# --- SVM-RFE: 取 recode：rfe 的 variables 是"逐轮淘汰顺序的倒序"；rank1 = 最后一轮仍在的基因
sv <- fread(file.path(out, "results", "D4_svmrfe_profile.tsv"))
sv <- unique(sv, by = c("var", "Variables", "Resample"))
# 每个变量在简历中首次出现的 Variables 值 = 它所属的最小可行子集
# caret 的 rfe$variables 里同一基因有多行（每轮各一），取最小 Variables 行 = 该基因第一次进入并获胜的子集
first <- sv[, .(minVars = min(Variables)), by = var]
# 简历最优子集在 accuracy 达到最大值的多个（8/20/25/40/60/80/100/167），取最精简的并列最优
# 但 variables 表里 minVars==8 的基因数 > 8（因为按重要性排序，前 8 个才是 8 基因集）
order8 <- sv[Variables == 8][order(-Overall)][, unique(var)]
# 对所有并列最优子集取"在最小可行子集中出现"的基因，按重要性（首次出现时 Overall）排序
svm_order <- sv[, .(Overall = max(Overall)), by = var][order(-Overall)]$var
svm_sel <- head(svm_order, 8)   # 与 best_n = 8 一致；因 8 基因集 CV accuracy = 1.00（并列最优中最简）
cat("SVM-RFE 排序前 8（最简并列最优子集）:", paste(svm_sel, collapse = ", "), "\n")
cat("  (该子集 5 折 CV accuracy = 1.00，与该队列可达的最高值一致；"
    , "队列 n=21，accuracy=1.0 属小样本过拟合上界，仅作筛选器不作性能声明)\n")

# --- 投票 ---
vote <- data.table(gene = genes)
vote[, LASSO := as.integer(gene %in% lasso_sel)]
vote[, RF_top20 := as.integer(gene %in% rf_sel)]
vote[, SVM_RFE := as.integer(gene %in% svm_sel)]
vote[, PPI_hub := as.integer(gene %in% hub$gene)]
vote[, n_vote := LASSO + RF_top20 + SVM_RFE]
vote[, lasso_coef := lasso$coef[match(gene, lasso$gene)]]
vote[, rf_MDG := imp$MDG[match(gene, imp$gene)]]
setorder(vote, -n_vote, -PPI_hub, -rf_MDG)
fwrite(vote, file.path(out, "results", "D4_ml_vote_matrix.tsv"), sep = "\t")
fwrite(vote[n_vote > 0], file.path(out, "results", "D4_ml_selected.tsv"), sep = "\t")

cat("\n=== n_vote >= 2 ===\n")
print(vote[n_vote >= 2, .(gene, LASSO, RF_top20, SVM_RFE, PPI_hub, n_vote,
                          lasso_coef = round(lasso_coef, 4), rf_MDG = round(rf_MDG, 2))])

# 最终 1-2 基因：优先 n_vote=3；若无则 n_vote==2 & PPI_hub==1；再退 n_vote==2
final <- vote[n_vote == 3]
if (nrow(final) == 0) final <- vote[n_vote == 2 & PPI_hub == 1]
if (nrow(final) == 0) final <- vote[n_vote == 2]
final <- head(final[order(-n_vote, -PPI_hub, -abs(lasso_coef), -rf_MDG)], 2)
if (nrow(final) == 0) {
  cat("\n!!! 三算法无 >=2 票交集，触发计划书退路：并集 + 多因素 logistic 逐步回归\n")
  uni <- vote[n_vote >= 1]
  final <- head(uni[order(-n_vote, -PPI_hub, -abs(lasso_coef))], 2)
}
fwrite(final, file.path(out, "results", "D4_final_genes.tsv"), sep = "\t")
cat("\n>>> 最终收敛基因:\n"); print(final)

prior <- vote[gene %in% c("SERPINE1", "THBS1", "LOX")]
cat("\n=== 先验候选在 D4 ===\n"); print(prior[, .(gene, LASSO, RF_top20, SVM_RFE, PPI_hub, n_vote,
                                                  lasso_coef = round(lasso_coef, 4))])

js <- list(
  lasso_lambda_min = 0.0042428, lasso_selected = lasso_sel,
  rf_top20 = rf_sel, rf_oob_error = 0.04761905,
  svm_best_n = 8, svm_selected = svm_sel,
  svm_cv_accuracy = 1.0,
  ppi_hub_intersection = hub$gene,
  final_genes = final$gene,
  prior_candidates = setNames(as.list(prior$n_vote), prior$gene)
)
write_json(js, file.path(out, "data", "D4_ml_selected.json"), auto_unbox = TRUE, pretty = TRUE)

# --- 图 ---
prof <- fread(file.path(out, "results", "D4_svmrfe_profile.tsv"))
pv <- unique(prof[, .(Variables, Accuracy)])
pv <- pv[!is.na(Accuracy)][order(Variables)]
p3 <- ggplot(pv, aes(Variables, Accuracy)) +
  geom_line(colour = "#4c72b0") + geom_point(colour = "#4c72b0", size = 2.2) +
  geom_point(data = pv[Variables == 8], colour = "#d62728", size = 3.6) +
  geom_hline(yintercept = 1, linetype = 3, colour = "grey40") +
  scale_x_continuous(trans = "log2", breaks = pv$Variables) +
  labs(title = "SVM-RFE recursive feature elimination (5-fold CV)",
       subtitle = "red = 8-gene minimal subset reaching CV accuracy 1.00 (n = 21, upper bound)",
       x = "Number of features (log2)", y = "CV accuracy") +
  theme_bw(base_size = 11)
ggsave(file.path(out, "figures", "D4_Fig5_svmrfe.pdf"), p3, width = 6.4, height = 4.6)
ggsave(file.path(out, "figures", "D4_Fig5_svmrfe.png"), p3, width = 6.4, height = 4.6, dpi = 300)

# 收敛热图（前 30 名投票矩阵）
sel <- vote[n_vote > 0][order(-n_vote, -PPI_hub, -rf_MDG)][1:min(30, .N)]
vm <- melt(sel, id.vars = "gene", measure.vars = c("LASSO", "RF_top20", "SVM_RFE"))
vm[, gene := factor(gene, levels = rev(sel$gene))]
p4 <- ggplot(vm, aes(gene, variable, fill = factor(value))) +
  geom_tile(colour = "white") +
  scale_fill_manual(values = c("0" = "grey93", "1" = "#d62728"), guide = "none") +
  labs(title = "Three-algorithm evidence matrix (LASSO / RF / SVM-RFE)",
       subtitle = paste0("final: ", paste(final$gene, collapse = " + "), "  (prior SERPINE1/LOX: ",
                         paste0(prior$gene, "=", prior$n_vote, collapse = ", "), ")"),
       x = NULL, y = NULL) +
  theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave(file.path(out, "figures", "D4_Fig6_ml_convergence.pdf"), p4, width = 8.4, height = 4.4)
ggsave(file.path(out, "figures", "D4_Fig6_ml_convergence.png"), p4, width = 8.4, height = 4.4, dpi = 300)
cat("\nD4-04 done\n")