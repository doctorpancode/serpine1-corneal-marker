# ---------------------------------------------------------------------------
# 排名稳定性：20 次独立网络重建 + 20 次度保持置换零分布
#   输入：第8天执行结果/data/D8_knk_qcCounts.rds （947 x 3587 CPM）
#   复现 03_regularized_knockout.R 的 mean|Δ| 与 loess(out_degree) 残差口径
#   输出：results/D9_rank_stability.tsv（重建排名）
#         results/D9_null_ranks.tsv   （置换零分布排名）
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(scTenifoldNet)
  library(scTenifoldKnk)
  library(igraph)
})

D8  <- "D:/A生信分析/虚拟敲除/第8天执行结果"
RES <- "D:/A生信分析/虚拟敲除/第9天执行结果/results"
CORE <- c("Serpine1","Tmsb4x","P4hb","S100a13","Sod1","Prdx2","Vdac1","Capg","Gja1","Xiap","Vim")
EPI  <- c("Krt12","Aqp5","Gsta4","Upk3bl","Upk1b","Krt15","Krt17","Krt19","Krt13",
          "Dsg1a","Perp","Tmprss11d","Lypd3","Calml3","Gsta2","Sfn","Krt5","Krt14")
NREP <- 20
HEAT_T <- 10
SPAN <- 0.5

qcX <- readRDS(file.path(D8, "data", "D8_knk_qcCounts.rds"))
cat("input:", dim(qcX), "\n")

score <- function(WT) {
  WT <- as.matrix(WT)
  diag(WT) <- 0
  genes <- rownames(WT)
  H <- scTenifoldKnk::heatKernel(WT, t = HEAT_T)
  S <- scTenifoldKnk:::.logCPM(as(qcX, "dgCMatrix"), genes)
  n <- length(genes)
  E <- matrix(NA_real_, n, n)
  for (i in seq_len(n)) {
    E[i, ] <- abs(as.numeric(crossprod(-S$mu[i] / S$sd[i], H[i, ])) * S$sd)
  }
  diag(E)[] <- NA
  dimnames(E) <- list(genes, genes)
  epi2 <- intersect(EPI, genes)
  dt <- data.table(ko_gene = genes, out_degree = as.numeric(rowSums(WT != 0)))
  dt[, mean_raw := rowMeans(E[, epi2], na.rm = TRUE)]
  dt[, resid := resid(loess(mean_raw ~ out_degree, data = dt, span = SPAN))]
  dt[, rank_raw := frank(-mean_raw, ties.method = "first")]
  dt[, rank_resid := frank(-resid, ties.method = "first")]
  dt
}

## ---------------- A. 20 次重建 ----------------
stab <- list()
for (s in seq_len(NREP)) {
  set.seed(20261110 + s)
  t0 <- Sys.time()
  nets <- scTenifoldNet::makeNetworks(qcX, nNet = 10, nCells = 500, nComp = 3,
                                      scaleScores = TRUE, symmetric = FALSE, q = 0.9)
  td <- scTenifoldNet::tensorDecomposition(nets, K = 3, maxError = 1e-5,
                                           maxIter = 1000, nDecimal = 3)
  WT <- td$X
  dt <- score(WT)
  dt[, rep := s]
  stab[[s]] <- dt[, .(rep, ko_gene, out_degree, mean_raw, resid, rank_raw, rank_resid)]
  cat(sprintf("rebuild %2d  %.1f s  Tmsb4x rank_resid=%d  Serpine1=%d\n", s,
              as.numeric(difftime(Sys.time(), t0, units = "secs")),
              dt[ko_gene == "Tmsb4x", rank_resid], dt[ko_gene == "Serpine1", rank_resid]))
}
STAB <- rbindlist(stab)
fwrite(STAB, file.path(RES, "D9_rank_stability.tsv"), sep = "\t")

sum_stab <- STAB[ko_gene %in% CORE, .(
  n_rep = .N,
  rank_raw_med = median(rank_raw), rank_raw_min = min(rank_raw), rank_raw_max = max(rank_raw),
  rank_resid_med = median(rank_resid), rank_resid_min = min(rank_resid),
  rank_resid_max = max(rank_resid),
  n_top1_resid = sum(rank_resid == 1), n_top10_resid = sum(rank_resid <= 10),
  n_bottom_third = sum(rank_resid > 631)
), by = ko_gene][order(rank_resid_med)]
fwrite(sum_stab, file.path(RES, "D9_rank_stability_summary.tsv"), sep = "\t")
print(sum_stab)

## ---------------- B. 20 次度保持置换零分布 ----------------
WT0 <- as.matrix(readRDS(file.path(D8, "data", "D8_WT_network.rds")))
diag(WT0) <- 0
g0 <- graph_from_adjacency_matrix(WT0 != 0, mode = "directed")
cat("edges:", ecount(g0), " density:", round(ecount(g0) / (vcount(g0) * (vcount(g0) - 1)), 4), "\n")

nulls <- list()
for (s in seq_len(NREP)) {
  set.seed(20261110 + 100 + s)
  t0 <- Sys.time()
  gp <- rewire(g0, keeping_degseq(niter = ecount(g0) * 2))
  Wp <- as.matrix(as_adjacency_matrix(gp)) * 1
  # 保留原网络的边权分布：按非零位置随机重排权重
  set.seed(20261110 + 200 + s)
  Wp[Wp != 0] <- sample(WT0[WT0 != 0], sum(Wp != 0), replace = TRUE)
  rownames(Wp) <- colnames(Wp) <- rownames(WT0)
  dt <- score(Wp)
  dt[, rep := s]
  nulls[[s]] <- dt[, .(rep, ko_gene, out_degree, mean_raw, resid, rank_raw, rank_resid)]
  cat(sprintf("perm  %2d  %.1f s  Tmsb4x rank_resid=%d  Serpine1=%d\n", s,
              as.numeric(difftime(Sys.time(), t0, units = "secs")),
              dt[ko_gene == "Tmsb4x", rank_resid], dt[ko_gene == "Serpine1", rank_resid]))
}
NULLS <- rbindlist(nulls)
fwrite(NULLS, file.path(RES, "D9_null_ranks.tsv"), sep = "\t")

sum_null <- NULLS[ko_gene %in% CORE, .(
  n_rep = .N,
  rank_resid_med = median(rank_resid), rank_resid_min = min(rank_resid),
  rank_resid_max = max(rank_resid), n_top1_resid = sum(rank_resid == 1)
), by = ko_gene][order(rank_resid_med)]
fwrite(sum_null, file.path(RES, "D9_null_ranks_summary.tsv"), sep = "\t")
print(sum_null)

cat("\nDONE\n")