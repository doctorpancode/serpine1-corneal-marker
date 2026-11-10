# ============================================================
# D7-02c 使用 scTenifoldKnk 2.0.0 官方函数重做虚拟敲除
#
# 流程完全按官方 pipeline：
#   1. scQC              质量控制
#   2. cpmNormalization  CPM 归一化
#   3. makeNetworks      nNet=10 次下采样构建 PC-regression 网络
#   4. tensorDecomposition 张量分解（K=3）合成单一 WT 网络
#   5. strictDirection   （包内部函数）
#   6. 敲除：WT[gKO, ] <- 0
#   7. manifoldAlignment 流形对齐
#   8. dRegulation       差异调控分析（含经验零分布）
#
# 关键约束（官方文档）：tensorDecomposition 是 dense O(n^2) 内存，
#   基因数必须压到 ~500 量级。我们用「D3 池 + 核心基因 + 高变基因」
#   优先级排序取前 N，确保 10 个核心基因与 Serpine1/Tmsb4x 全在网内。
# ============================================================
suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(Matrix)
  library(scTenifoldKnk); library(scTenifoldNet)
})
set.seed(20261109)
out <- "D:/A生信分析/虚拟敲除/第7天执行结果"
d6  <- "D:/A生信分析/虚拟敲除/第6天执行结果"

N_GENES <- 500      # 官方推荐量级
run_tag <- "D7_knk_official"

sc <- readRDS(file.path(d6, "data", "D6_seurat_GSE275678_final.rds"))
sc <- JoinLayers(sc)
cat("cells:", ncol(sc), " genes:", nrow(sc), "\n")

# ---- 1. 选基因：核心基因 > D3 池 > 高变基因 ----
pool <- fread("D:/A生信分析/虚拟敲除/第3天执行结果/data/candidate_pool.csv")
p2m <- function(x) paste0(substr(x, 1, 1), tolower(substr(x, 2, nchar(x))))
core <- c("Serpine1", "Tmsb4x", "P4hb", "S100a13", "Sod1", "Prdx2", "Vdac1", "Capg", "Gja1", "Xiap")
hv <- VariableFeatures(sc)
pool_m <- intersect(unique(p2m(pool$gene)), rownames(sc))

genes <- unique(c(intersect(core, rownames(sc)), pool_m, head(hv, 2000)))
# 表达过滤：至少 5% 细胞检出（scQC 默认 minPCT=0.05）
cnt_all <- GetAssayData(sc, layer = "counts")
pct <- Matrix::rowSums(cnt_all[, 1:4000] > 0) / 4000
genes <- genes[genes %in% names(pct)[pct >= 0.05]]
genes <- head(genes, N_GENES)
# 保证核心基因在首位
genes <- unique(c(intersect(core, genes), genes))
genes <- genes[genes %in% rownames(cnt_all)]
cat("genes selected:", length(genes), "\n")

# ---- 2. 观测矩阵（官方要求 raw counts, genes x cells）----
X <- as.matrix(cnt_all[genes, ])
cat("matrix:", nrow(X), "genes x", ncol(X), "cells; nnz =",
    round(100 * sum(X > 0) / length(X), 1), "%\n")
rm(cnt_all); gc()

# 下采样细胞到 5000（保持 4 位小数内存可控，nCells 参数会再抽样）
if (ncol(X) > 5000) X <- X[, sample(ncol(X), 5000)]
cat("cells used:", ncol(X), "\n")
saveRDS(X, file.path(out, "data", paste0(run_tag, "_input.rds")))

# ---- 3. 跑官方 pipeline ----
t0 <- Sys.time()
res <- scTenifoldKnk(
  countMatrix   = X,
  gKO           = c("Serpine1", "Tmsb4x"),
  transcriptomeWide = FALSE,
  qc            = TRUE,
  qc_minLibSize = 1000,
  qc_removeOutlierCells = TRUE,
  qc_minPCT     = 0.05,
  qc_maxMTratio = 0.1,
  nc_nNet       = 10,
  nc_nCells     = 500,
  nc_nComp      = 3,
  nc_scaleScores = TRUE,
  nc_q          = 0.9,
  td_K          = 3,
  td_maxIter    = 1000,
  td_maxError   = 1e-5,
  td_nDecimal   = 3,
  ma_nDim       = 2,
  ma_method     = "manifold",
  ma_heatT      = 10,
  dr_empiricalNull = TRUE,
  dr_direction  = TRUE,
  dr_directionT = 5,
  nCores        = 4,
  seed          = 20261109
)
cat("pipeline done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

# ---- 4. 落盘 ----
saveRDS(res, file.path(out, "data", paste0(run_tag, "_result.rds")))

wt <- res$tensorNetworks$WT
cat("WT network dim:", nrow(wt), "x", ncol(wt), "; nonzero edges:", sum(wt != 0), "\n")
for (g in c("Serpine1", "Tmsb4x")) {
  if (g %in% rownames(wt)) cat(sprintf("  %s outdegree nonzero = %d\n", g, sum(wt[g, ] != 0)))
}

dr <- as.data.table(res$diffRegulation)
fwrite(dr, file.path(out, "results", paste0(run_tag, "_diffRegulation.tsv")), sep = "\t")
cat("diffRegulation rows:", nrow(dr), " cols:", paste(names(dr), collapse = ", "), "\n")

# 方向表
if (!is.null(res$perturbationDirections)) {
  print(dim(res$perturbationDirections))
} else if (!is.null(res$manifoldAlignment)) {
  cat("manifoldAlignment dim:", dim(res$manifoldAlignment), "\n")
}

cat("\n=== top differentially regulated (all rows) ===\n")
print(head(dr[order(p.value)], 15))
cat("\nD7-02c done\n")