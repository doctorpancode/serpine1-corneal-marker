# ============================================================
# D7-02 scTenifoldKnk 虚拟敲除：Serpine1 / Tmsb4x
# 输入：D6 主队列 final 对象（GSE275678，11 细胞类型）
# 策略：取全部细胞构建 PC-regression 网络（与 scTenifoldNet 同源算法），
#       对两个主干基因分别做 in-silico knockout，比较扰动基因集
# ============================================================
suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(Matrix)
})
set.seed(20261109)
out <- "D:/A生信分析/虚拟敲除/第7天执行结果"
d6  <- "D:/A生信分析/虚拟敲除/第6天执行结果"

sc <- readRDS(file.path(d6, "data", "D6_seurat_GSE275678_final.rds"))
sc <- JoinLayers(sc)
cat("cells:", ncol(sc), " genes:", nrow(sc), "\n")

# ---- 1. 取高变基因（scTenifoldKnk 建议 1000-3000 基因）+ 强制纳入两个主干与池内基因 ----
pool <- fread("D:/A生信分析/虚拟敲除/第3天执行结果/data/candidate_pool.csv")
core <- c("Serpine1", "Tmsb4x", "P4hb", "S100a13", "Sod1", "Prdx2", "Vdac1", "Capg", "Gja1", "Xiap")
hv <- VariableFeatures(sc)
p2m <- function(x) paste0(substr(x, 1, 1), tolower(substr(x, 2, nchar(x))))
pool_m <- intersect(unique(p2m(pool$gene)), rownames(sc))
genes <- unique(c(head(hv, 2000), intersect(core, rownames(sc)), pool_m))
genes <- intersect(genes, rownames(sc))
cat("network genes:", length(genes), "\n")

X <- as.matrix(GetAssayData(sc, layer = "data")[genes, ])
X <- X[rowSums(X > 0) >= 0.03 * ncol(X), , drop = FALSE]
cat("after expression filter:", nrow(X), "genes x", ncol(X), "cells\n")

# ---- 2. 细胞下采样（PC regression 内存 O(n*p)）----
if (ncol(X) > 6000) {
  idx <- sample(ncol(X), 6000)
  X <- X[, idx]
}
cat("cells used:", ncol(X), "\n")

# ---- 3. PC regression 邻接矩阵（scTenifoldNet 的核心）----
t0 <- Sys.time()
pcNet <- function(x, nComp = 3, nCores = 1) {
  # x: genes x cells
  y <- t(x)                                  # cells x genes
  y <- scale(y)
  pr <- prcomp(y, center = FALSE, scale. = FALSE, rank. = nComp)
  pcs <- pr$x[, 1:nComp, drop = FALSE]       # cells x nComp
  W <- sapply(seq_len(nrow(x)), function(i) {
    fit <- lm.fit(x = pcs, y = y[, i])
    coefs <- fit$coefficients
    if (any(is.na(coefs))) coefs[is.na(coefs)] <- 0
    coefs
  })                                          # nComp x genes
  S <- solve(crossprod(pcs) / ncol(y) + diag(1e-6, nComp))
  tW <- t(W)
  A <- tW %*% S %*% W     # genes x genes
  A
}
W <- pcNet(X, nComp = 3)
cat("network built in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "min\n")
rownames(W) <- colnames(W) <- rownames(X)
cat("edges (nonzero):", sum(W != 0), "\n")
saveRDS(W, file.path(out, "data", "D7_pcreg_network.rds"))

# 对称化 + 归一
A <- (W + t(W)) / 2
diag(A) <- 0
saveRDS(A, file.path(out, "data", "D7_network_sym.rds"))

# ---- 4. in-silico knockout：把目标基因的边去掉 ----
ko <- function(A, gene) {
  if (!gene %in% rownames(A)) return(NULL)
  B <- A
  B[gene, ] <- 0; B[, gene] <- 0
  B
}

run_ko <- function(target, nPerm = 20) {
  B <- ko(A, target)
  if (is.null(B)) return(NULL)
  # 拓扑差异：每个基因的度与特征向量中心性变化
  d_orig <- Matrix::rowSums(A != 0)
  d_ko   <- Matrix::rowSums(B != 0)
  c_orig <- tryCatch(as.numeric(igraph::evcent(igraph::graph_from_adjacency_matrix(A, weighted = TRUE))$vector), error = function(e) rep(NA, nrow(A)))
  c_ko   <- tryCatch(as.numeric(igraph::evcent(igraph::graph_from_adjacency_matrix(B, weighted = TRUE))$vector), error = function(e) rep(NA, nrow(A)))
  names(d_orig) <- names(d_ko) <- rownames(A)
  names(c_orig) <- names(c_ko) <- rownames(A)
  delta_deg <- d_orig - d_ko
  delta_cen <- c_orig - c_ko
  d0 <- data.table(gene = rownames(A), target = target,
                   degree = as.numeric(d_orig), delta_degree = as.numeric(delta_deg),
                   centrality = c_orig, delta_centrality = c_orig - c_ko)
  # 置换检验：随机敲除同度数的基因，得到 delta_degree 的零分布
  deg_rank <- rank(-d_orig)
  cand <- names(d_orig)[d_orig > 0]
  null <- replicate(nPerm, {
    set.seed(NULL)
    g <- sample(cand, 1)
    Bn <- A; Bn[g, ] <- 0; Bn[, g] <- 0
    d_next <- Matrix::rowSums(Bn != 0)
    sum(d_orig - d_next)
  })
  d0[, `:=`(target_degree = as.numeric(d_orig[target]),
            perm_mean = mean(null), perm_max = max(null))]
  d0
}

res <- rbindlist(lapply(c("Serpine1", "Tmsb4x"), run_ko))
fwrite(res, file.path(out, "results", "D7_knockout_delta.tsv"), sep = "\t")

cat("\nSerpine1 KO:\n")
print(res[target == "Serpine1"][order(-delta_degree)][1:15, .(gene, degree, delta_degree = round(delta_degree, 2))])
cat("\nTmsb4x KO:\n")
print(res[target == "Tmsb4x"][order(-delta_degree)][1:15, .(gene, degree, delta_degree = round(delta_degree, 2))])

# ---- 5. 扰动基因集与 D3 候选池交集 ----
pert <- res[delta_degree != 0]
pert[, in_pool := toupper(gene) %in% pool$gene]
fwrite(pert, file.path(out, "results", "D7_knockout_perturbed.tsv"), sep = "\t")
for (t in c("Serpine1", "Tmsb4x")) {
  pp <- pert[target == t]
  cat("\n", t, ": perturbed", nrow(pp), "genes; in D3 pool",
      sum(pp$in_pool), "(", round(100 * mean(pp$in_pool), 1), "%)\n")
}

cat("\nD7-02 done\n")