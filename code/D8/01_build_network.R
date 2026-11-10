# ============================================================
# D8-01 扩网络：把基因数从 D7 的 494 提到官方上限量级(~900)
#
# 背景（D7 的两条硬约束，都来自 scTenifoldKnk 2.0.0 源码，不是经验）：
#   1) scQC() 里的基因过滤是 `rowMeans(X != 0) > minPCT`，对全细胞算，
#      不是 D7-02c 里那种 `[, 1:4000]` 近似。D7-02c 用前 4000 个细胞粗筛
#      (>=5%)，再让官方 scQC 严格筛，500 -> 494 掉 6 个。
#      本脚本直接先用官方函数 scQC() 做一次预筛，选基因时保证 100% 存活。
#   2) tensorDecomposition 是 dense O(n^2) x nNet=10 张 494x494 的矩阵，
#      内存约 10 * 8 * n^2 字节 + CP 分解工作区。
#
# 所以 N_GENES 先按 checkMemory() 的实测内存模型算出可用上限，再取
# min(目标值, 上限 * 0.85) 保险。
#
# 输出：
#   data/D8_knk_input.rds            原始 counts（genes x cells，官方要求）
#   data/D8_knk_qcCounts.rds         scQC + CPM 之后的矩阵（03 复用它，保证同网络）
#   data/D8_WT_network.rds           张量分解后的 WT 邻接矩阵（900x900，511 MB 以内）
#   results/D8_network_gene_list.tsv 入选基因 + 来源标记
#   results/D8_network_edge_summary.tsv 每个核心基因的出/入度
# ============================================================
suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(Matrix)
  library(scTenifoldKnk); library(scTenifoldNet)
})
set.seed(20261110)

out <- "D:/A生信分析/虚拟敲除/第8天执行结果"
d6  <- "D:/A生信分析/虚拟敲除/第6天执行结果"
d3  <- "D:/A生信分析/虚拟敲除/第3天执行结果"
dir.create(file.path(out, "data"),    showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out, "results"), showWarnings = FALSE, recursive = TRUE)

N_TARGET <- 950          # 目标基因数（官方 vignette 上限 1000，留 5% 冗余）
N_CELLS  <- 5000         # 与 D7 一致，保证两次 WT 网络之间只差基因数
CORE     <- c("Serpine1", "Tmsb4x", "P4hb", "S100a13", "Sod1", "Prdx2",
              "Vdac1", "Capg", "Gja1", "Xiap", "Vim", "Cdh2")

# ---------------- 1. 读单细胞主队列 ----------------
sc <- readRDS(file.path(d6, "data", "D6_seurat_GSE275678_final.rds"))
sc <- JoinLayers(sc)
cnt <- GetAssayData(sc, layer = "counts")
hv  <- VariableFeatures(sc)          # 必须在 rm(sc) 前取，否则变量丢失
cat("GSE275678:", nrow(cnt), "genes x", ncol(cnt), "cells\n")
rm(sc); gc()

# 下采样细胞（再抽 8000 给 scQC 里 outlier 剔除留余量）
if (ncol(cnt) > 8000) cnt <- cnt[, sample(ncol(cnt), 8000)]

# ---------------- 2. 官方式预筛：先跑一次 scQC 看谁真能活下来 ----------------
# 坑：scTenifoldKnk::scQC 的 MT 过滤写法是 X[mtGenes, ]，在 dgCMatrix 上
#     会因 drop=TRUE 退化成数值向量，随后 colSums 报
#     "'x' must be an array of at least two dimensions"。
#     D7 的 500 基因恰好不含 MT- 开头的基因，所以没踩到；D8 的 950 个里有 1 个。
#     这里直接用 Matrix 版的等价实现，过滤条件与官方函数逐条一致。
scQC_safe <- function(X, minLibSize = 1000, removeOutlierCells = TRUE,
                      minPCT = 0.05, maxMTratio = 0.1) {
  X <- as(X, "dgCMatrix")
  X@x[X@x < 0] <- 0
  X <- X[, Matrix::colSums(X) > minLibSize]
  if (removeOutlierCells) {
    ls <- Matrix::colSums(X)
    X <- X[, !ls %in% boxplot.stats(ls)$out]
  }
  mt <- grepl("^MT-", toupper(rownames(X)), ignore.case = TRUE)
  if (any(mt)) {
    mtr <- Matrix::colSums(X[mt, , drop = FALSE]) / Matrix::colSums(X)
    X <- X[, mtr < maxMTratio]
  }
  X <- X[Matrix::rowMeans(X != 0) > minPCT, ]
  as(X, "dgCMatrix")
}
qc0 <- scQC_safe(cnt)
survivors <- rownames(qc0)
cat("survivors after scQC:", length(survivors), "x", ncol(qc0), "\n")
rm(qc0); gc()

# ---------------- 3. 基因优先级 ----------------
pool <- fread(file.path(d3, "data", "candidate_pool.csv"))
p2m <- function(x) paste0(substr(x, 1, 1), tolower(substr(x, 2, nchar(x))))
pool[, mouse := p2m(gene)]
pool[, tier_rank := fifelse(tier == "T1_triple_DEG", 1L,
                     fifelse(tier == "T2_double", 2L, 3L))]
setorder(pool, tier_rank, -n_source_sets)
pool_m <- unique(pool$mouse)

# 核心基因优先保留（含 D8 要敲的 P4hb / Vim / S100a13）
core_ok <- intersect(CORE, survivors)
cat("core genes present:", paste(core_ok, collapse = ", "), "\n")

hv <- head(hv[hv %in% survivors], 3000)

rank_gene <- function(g) {
  # 1=核心, 2=D3池, 3=高变, 4=按平均表达
  if (g %in% core_ok) return(1L)
  if (g %in% pool_m) return(2L)
  if (g %in% hv) return(3L)
  4L
}
cand <- survivors
rk <- vapply(cand, rank_gene, integer(1))
# 同一优先级内：D3 池按 tier/n_source_sets，其余按平均表达量降序
detv <- as.numeric(Matrix::rowMeans(cnt[cand, , drop = FALSE] > 0))
pri <- data.table(gene = cand, grp = rk, det = detv)
pri[, in_pool := gene %in% pool_m]
pri <- merge(pri, pool[, .(gene = mouse, tier_rank, n_source_sets)],
             by = "gene", all.x = TRUE, sort = FALSE)
pri[is.na(tier_rank), tier_rank := 99L]
pri[is.na(n_source_sets), n_source_sets := 0L]
setorder(pri, grp, tier_rank, -n_source_sets, -det)
genes <- pri$gene

# ---------------- 4. 内存上限裁决 ----------------
mem <- scTenifoldNet:::checkMemory(length(genes), nNet = 10, nConditions = 1, warn = FALSE)
cat(sprintf("memory: available %.1f GB / total %.1f GB; required for %d genes = %.1f GB; maxGenes = %s\n",
            mem$available/1e9, mem$total/1e9, length(genes), mem$required/1e9,
            ifelse(is.na(mem$maxGenes), "NA", mem$maxGenes)))
cap_genes <- if (is.na(mem$maxGenes)) N_TARGET else min(N_TARGET, floor(mem$maxGenes * 0.85))
cat("gene cap from memory:", cap_genes, "\n")

# 核心基因必须在内：先放核心，再按优先级补齐
keep <- unique(c(intersect(core_ok, genes), head(genes, cap_genes)))
keep <- keep[keep %in% survivors]
cat("genes selected:", length(keep), "\n")
stopifnot(all(core_ok %in% keep))

fwrite(pri[gene %in% keep], file.path(out, "results", "D8_network_gene_list.tsv"), sep = "\t")

# ---------------- 5. 构造矩阵并落盘 ----------------
X <- as.matrix(cnt[keep, ])
if (ncol(X) > N_CELLS) X <- X[, sample(ncol(X), N_CELLS)]
cat("matrix:", nrow(X), "genes x", ncol(X), "cells; nnz =",
    round(100 * sum(X > 0)/length(X), 1), "%\n")
saveRDS(X, file.path(out, "data", "D8_knk_input.rds"))

# ---------------- 6. 与 D7 完全一致地重建 WT 网络 ----------------
t0 <- Sys.time()
qcX <- scQC_safe(X)
cat("after scQC:", nrow(qcX), "x", ncol(qcX), "\n")
qcX <- as(qcX, "dgCMatrix")
cpmX <- scTenifoldNet::cpmNormalization(qcX)
saveRDS(cpmX, file.path(out, "data", "D8_knk_qcCounts.rds"))

set.seed(20261110)
netList <- scTenifoldNet::makeNetworks(X = cpmX, q = 0.9, nNet = 10, nCells = 500,
                                       scaleScores = TRUE, symmetric = FALSE,
                                       nComp = 3, nCores = 6)
cat("networks built in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

set.seed(20261110)
td <- scTenifoldNet::tensorDecomposition(xList = netList, K = 3, maxError = 1e-5,
                                         maxIter = 1000, nDecimal = 3, seed = 20261110)
WT <- td$X
WT <- scTenifoldKnk:::strictDirection(WT, lambda = 0)
WT <- as.matrix(WT)
diag(WT) <- 0
WT <- t(WT)                                  # 官方：转置后行为"出边"
rm(netList, td); gc()

cat("WT network:", nrow(WT), "x", ncol(WT), "; nonzero edges:", sum(WT != 0),
    "; density =", round(100 * sum(WT != 0)/nrow(WT)^2, 2), "%\n")
saveRDS(WT, file.path(out, "data", "D8_WT_network.rds"))

# ---------------- 7. 出/入度体检 ----------------
cnt_genes <- CORE
inn <- cnt_genes %in% rownames(WT)
g2 <- cnt_genes[inn]
es <- data.table(
  gene       = cnt_genes,
  in_network = inn,
  out_degree = NA_integer_,
  in_degree  = NA_integer_
)
es[inn, out_degree := as.integer(rowSums(WT[g2, , drop = FALSE] != 0))]
es[inn, in_degree  := as.integer(colSums(WT[, g2, drop = FALSE] != 0))]
print(es)
fwrite(es, file.path(out, "results", "D8_network_edge_summary.tsv"), sep = "\t")

# ---------------- 8. 元数据 ----------------
meta <- data.table(
  item = c("input_cells_total", "cells_used_for_networks", "scQC_genes", "scQC_cells",
           "N_TARGET", "maxGenes_by_memory", "norm_explained_note",
           "WT_genes", "WT_nonzero_edges", "WT_density_pct", "missing_core"),
  value = c(15873, N_CELLS, nrow(qcX), ncol(qcX),
            N_TARGET, as.character(mem$maxGenes),
            "见 logs/01_build_network.txt 的 CP decomposition 行",
            nrow(WT), sum(WT != 0), round(100 * sum(WT != 0)/nrow(WT)^2, 2),
            paste(setdiff(CORE, rownames(WT)), collapse = ","))
)
fwrite(meta, file.path(out, "metadata", "D8_network_meta.tsv"), sep = "\t")

cat(sprintf("\ntotal %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat("D8-01 done\n")