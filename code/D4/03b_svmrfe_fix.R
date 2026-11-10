# ============================================================
# D4-03b SVM-RFE 修复版（caret + svmRadial，直接用 caret 的 ROC 指标名）
# 单独运行，避免重跑整个 03 脚本
# ============================================================
suppressPackageStartupMessages({
  library(Biobase); library(data.table); library(caret)
})
set.seed(20261109)
out <- "D:/A生信分析/虚拟敲除/第4天执行结果"
es <- readRDS("D:/A生信分析/虚拟敲除/D1执行结果/D1_CNV_角膜新生血管_数据分析/data/processed/GSE241715_expression.rds")
x <- exprs(es); pd <- pData(es)
pool <- fread("D:/A生信分析/虚拟敲除/第3天执行结果/data/candidate_pool.csv")
genes <- intersect(pool$gene, rownames(x))
y <- factor(ifelse(pd$group == "Control", "Normal", "Lesion"), levels = c("Normal", "Lesion"))
X <- t(x[genes, , drop = FALSE])
v <- apply(X, 2, var); X <- X[, v > 1e-8, drop = FALSE]
cs <- scale(X); cs[is.na(cs)] <- 0

# 方差不齐 → 用 preProcess 内的 center/scale 已够，svmRadial 默认 tuneLength=3 sigma=NA -> 用 CARE 网格
subsets <- c(1, 2, 3, 5, 8, 12, 16, 20, 25, 30, 40, 50, 60, 80, 100, ncol(cs))
subsets <- sort(unique(subsets[subsets <= ncol(cs)]))

ctrl <- rfeControl(functions = caretFuncs, method = "cv", number = 5,
                   verbose = FALSE, allowParallel = FALSE)
set.seed(20261109)
rfe_fit <- rfe(x = as.data.frame(cs), y = y, sizes = subsets,
               rfeControl = ctrl, method = "svmRadial",
               preProcess = c("center", "scale"),
               trControl = trainControl(method = "cv", number = 5, classProbs = TRUE,
                                        summaryFunction = twoClassSummary, savePredictions = "final"),
               tuneLength = 4, metric = "ROC")
cat("best subset =", rfe_fit$bestSubset, "\n")
print(rfe_fit$results)
sv <- rfe_fit$variables
fwrite(as.data.table(sv), file.path(out, "results", "D4_svmrfe_profile.tsv"), sep = "\t")
best_n <- rfe_fit$bestSubset
svm_sel <- as.character(sv$var[sv$Variables == best_n])
cat("SVM-RFE best", best_n, "genes:", paste(svm_sel, collapse = ", "), "\n")
saveRDS(list(best_n = best_n, sel = svm_sel, results = rfe_fit$results),
        file.path(out, "data", "D4_svmrfe_fit.rds"))