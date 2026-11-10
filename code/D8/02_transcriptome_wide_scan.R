# ============================================================
# D8-02 全转录组规模的"敲除影响谱"扫描（热核批量对齐）
#
# 官方接口：scTenifoldKnk(transcriptomeWide = TRUE, ma_method = "heat")
#   -> hkManifoldAlignment(WT, qcCounts, gKO = as.list(targetGenes), t = ma_heatT)
#
# 为什么用热核而不是逐个 manifoldAlignment：
#   heatKernel(WT) 只算一次 900x900 的特征分解，之后每个基因的敲除效应
#   只是 H[g, ] 的一行线性组合（见 .kernelEffect），因此可以一次读进
#   几百个敲除，GPU/内存都省；manifoldAlignment 那套是每个基因重算一次
#   2 维嵌入，O(n) 次特征分解。
#
# 输入：
#   data/D8_WT_network.rds     01 产出的 900x900 WT 有向邻接
#   data/D8_knk_qcCounts.rds   01 产出的 scQC+CPM 矩阵（与建网时完全一致）
#
# 输出：
#   data/D8_scan_raw.rds       官方返回值（perturbationDistances / Directions）
#   results/D8_knockdown_spectrum.tsv       长表：ko_gene, affected_gene, distance, direction
#   results/D8_perturbation_distance_matrix.tsv  宽表：行=被敲基因，列=被扰动基因
#   results/D8_ko_impact_rank.tsv           每个敲除的整体影响强度排名
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(scTenifoldKnk); library(scTenifoldNet)
})
set.seed(20261110)
out <- "D:/A生信分析/虚拟敲除/第8天执行结果"

WT  <- readRDS(file.path(out, "data", "D8_WT_network.rds"))
qcX <- readRDS(file.path(out, "data", "D8_knk_qcCounts.rds"))
stopifnot(identical(rownames(WT), rownames(qcX)))
cat("WT:", nrow(WT), "x", ncol(WT), "; qcCounts:", nrow(qcX), "x", ncol(qcX), "\n")

HEAT_T <- 10                 # 官方 ma_heatT 默认值，与 D7-02c 一致
CORE   <- c("Serpine1", "Tmsb4x", "P4hb", "S100a13", "Sod1", "Prdx2",
            "Vdac1", "Capg", "Gja1", "Xiap", "Vim")
core_ok <- intersect(CORE, rownames(WT))

# ---------- 1. 热核批量对齐：全部基因逐个敲 ----------
# 官方 transcriptomeWide 分支内部就是 hkManifoldAlignment(WT, qcCounts, gKO, t)
# 这里直接调同名导出函数，跳过重复的建网步骤，同时保持与官方实现同一函数。
# 坑：hkManifoldAlignment 不接 H 之外的缓存时内部会自己算 heatKernel；
#     900x900 的特征分解约几秒，直接让它自己算，避免参数传递错位。
t0 <- Sys.time()
targets <- rownames(WT)                                  # 全 900+ 个基因逐个敲
D <- scTenifoldKnk::hkManifoldAlignment(WT, qcX, gKO = as.list(targets), t = HEAT_T)
cat("perturbation distance matrix:", nrow(D), "KO x", ncol(D), "genes\n")
cat("elapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "min\n")

set.seed(20261110)
S <- scTenifoldKnk::knockoutDirection(qcX, gKO = as.list(targets),
                                      genes = rownames(WT), t = 5)
S <- sign(as.matrix(S))
cat("direction matrix:", nrow(S), "x", ncol(S), "\n")

saveRDS(list(perturbationDistances = D, perturbationDirections = S,
             heatT = HEAT_T, nGenes = nrow(WT)),
        file.path(out, "data", "D8_scan_raw.rds"))

# ---------- 2. 落表 ----------
Dm <- as.data.table(as.table(D))
setnames(Dm, c("ko_gene", "affected_gene", "distance"))
Dm[, distance := as.numeric(distance)]
Dm[, ko_gene := as.character(ko_gene)]
Dm[, affected_gene := as.character(affected_gene)]
Sm <- as.data.table(as.table(S))
setnames(Sm, c("ko_gene", "affected_gene", "direction"))
Sm[, direction := as.numeric(direction)]
Sm[, ko_gene := as.character(ko_gene)]
Sm[, affected_gene := as.character(affected_gene)]

spec <- merge(Dm, Sm, by = c("ko_gene", "affected_gene"), all.x = TRUE)
spec[, abs_distance := abs(distance)]
spec[, self := ko_gene == affected_gene]
# 官方 hkManifoldAlignment 会把对角（自身）填成 logCPM 均值，这里单独标出，
# 下游算排名/相似度时剔除，避免"目标基因自己距离最大"把谱压平。
spec[self == TRUE, distance := NA_real_]
spec[self == TRUE, abs_distance := NA_real_]

setorder(spec, ko_gene, -abs_distance)
fwrite(spec[!is.na(distance)],
       file.path(out, "results", "D8_knockdown_spectrum.tsv"), sep = "\t")

# 宽表（行列同名，直接喂 ComplexHeatmap）
Dw <- D
Dw[cbind(seq_len(nrow(Dw)), match(rownames(Dw), colnames(Dw)))] <- NA
fwrite(data.table(ko_gene = rownames(Dw), as.data.table(Dw)),
       file.path(out, "results", "D8_perturbation_distance_matrix.tsv"), sep = "\t")

# ---------- 3. 每个敲除的影响强度排名 ----------
imp <- spec[!is.na(distance), .(
  n_genes      = .N,
  mean_abs     = mean(abs_distance),
  median_abs   = median(abs_distance),
  p90_abs      = quantile(abs_distance, 0.9),
  max_abs      = max(abs_distance),
  n_up         = sum(direction > 0, na.rm = TRUE),
  n_down       = sum(direction < 0, na.rm = TRUE),
  n_top10_hit  = NA_integer_
), by = ko_gene]
# top10 = 每个敲除自身谱里距离最大的 10 个基因，看它们是不是反复落在同一批基因上
top10 <- spec[!is.na(distance)][order(ko_gene, -abs_distance)][, head(.SD, 10), by = ko_gene]
shared <- top10[, .N, by = affected_gene][order(-N)]
cat("\n最常进入 top10 的受影响基因（敲除影响谱的公共靶点）：\n")
print(head(shared, 20))
fwrite(shared, file.path(out, "results", "D8_top10_shared_targets.tsv"), sep = "\t")
fwrite(imp[order(-mean_abs)], file.path(out, "results", "D8_ko_impact_rank.tsv"), sep = "\t")
fwrite(top10, file.path(out, "results", "D8_ko_top10_targets.tsv"), sep = "\t")

cat("\n=== 核心基因的谱宽度（mean |distance|） ===\n")
print(imp[ko_gene %in% core_ok][order(-mean_abs)])

# ---------- 4. 与 D7 的两敲结果对齐体检 ----------
# D7 用 494 基因子网 + manifoldAlignment，D8 用 900 基因网 + heatKernel。
# 方法不同不能逐值比较，但"受影响基因的排序"应当同向。这里只报 Spearman。
d7 <- fread("D:/A生信分析/虚拟敲除/第7天执行结果/results/D7_knk_official_twoKO.tsv")
chk <- list()
for (g in intersect(c("Serpine1", "Tmsb4x"), rownames(WT))) {
  a <- d7[target == g, .(gene, d7_dist = distance)]
  b <- spec[ko_gene == g, .(gene = affected_gene, d8_abs = abs_distance)]
  m <- merge(a, b, by = "gene")
  ct <- suppressWarnings(cor.test(m$d7_dist, m$d8_abs, method = "spearman"))
  chk[[g]] <- data.table(ko = g, n = nrow(m), rho = unname(ct$estimate), p = ct$p.value)
  # 各自 top10 的重叠
  t7 <- head(m[order(-d7_dist)]$gene, 10); t8 <- head(m[order(-d8_abs)]$gene, 10)
  chk[[g]][, overlap_top10 := paste(intersect(t7, t8), collapse = ",")]
}
chk <- rbindlist(chk)
print(chk)
fwrite(chk, file.path(out, "results", "D8_d7_d8_consistency.tsv"), sep = "\t")

cat("\nD8-02 done\n")