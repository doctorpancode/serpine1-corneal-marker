# ============================================================
# D7-02d 官方 scTenifoldKnk 逐基因敲除（替代一次跑两个基因）
# 理由：scTenifoldKnk(transcriptomeWide=FALSE) 一次只模拟一个 gKO 的
#       敲除；传两个基因等于同时敲除两个，diffRegulation 只给一份
#       合并结果，无法分别归因。
# 做法：固定同一 WT 网络（从 02c 的 rds 读回 tensorNetworks$WT），
#       分别对 Serpine1 / Tmsb4x 调 knockoutDirection + manifoldAlignment
#       + dRegulation（全部官方函数）。
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(scTenifoldKnk)
})
out <- "D:/A生信分析/虚拟敲除/第7天执行结果"
res <- readRDS(file.path(out, "data", "D7_knk_official_result.rds"))
WT  <- as.matrix(res$tensorNetworks$WT)
qcX <- readRDS(file.path(out, "data", "D7_knk_official_input.rds"))
cat("WT:", nrow(WT), "x", ncol(WT), "\n")

# 用与 02c 相同的 QC + CPM 复原 qcCounts（knockoutDirection 需要）
suppressPackageStartupMessages(library(scTenifoldNet))
qcX <- scQC(qcX, minLibSize = 1000, removeOutlierCells = TRUE,
            minPCT = 0.05, maxMTratio = 0.1)
qcX <- cpmNormalization(qcX)
cat("qcCounts:", nrow(qcX), "x", ncol(qcX), "\n")

allres <- list()
for (g in c("Serpine1", "Tmsb4x")) {
  cat("\n===== knocking out", g, "=====\n")
  if (!g %in% rownames(WT)) { cat("not in network\n"); next }
  KO <- WT
  KO[g, ] <- 0                      # 官方做法：只清出边
  cat("  WT outdegree:", sum(WT[g, ] != 0), "-> KO outdegree:", sum(KO[g, ] != 0), "\n")

  set.seed(20261109)
  dir <- tryCatch(
    scTenifoldKnk::knockoutDirection(qcX, gKO = list(g), genes = rownames(WT), t = 5),
    error = function(e) { cat("  knockoutDirection failed:", conditionMessage(e), "\n"); NULL })
  if (!is.null(dir)) dir <- setNames(as.numeric(dir), colnames(dir))

  set.seed(20261109)
  MA <- scTenifoldNet::manifoldAlignment(WT, KO, d = 2, nCores = 4)
  DR <- scTenifoldKnk::dRegulation(manifoldOutput = MA, gKO = g,
                                   empiricalNull = TRUE, direction = dir)
  dr <- as.data.table(DR)
  dr[, target := g]
  cat("  significant (FDR<0.05):", sum(dr$p.adj < 0.05, na.rm = TRUE), "\n")
  print(head(dr[order(p.value), .(gene, distance, Z, FC, p.value, p.adj, direction)], 15))
  allres[[g]] <- dr
}

ko2 <- rbindlist(allres, fill = TRUE)
fwrite(ko2, file.path(out, "results", "D7_knk_official_twoKO.tsv"), sep = "\t")
cat("\n===== 交集 =====")
s1 <- ko2[target == "Serpine1" & p.adj < 0.05]$gene
s2 <- ko2[target == "Tmsb4x"  & p.adj < 0.05]$gene
cat("\nSerpine1 KO sig (", length(s1), "):", paste(s1, collapse = ", "), "\n")
cat("Tmsb4x KO sig (", length(s2), "):", paste(s2, collapse = ", "), "\n")
cat("overlap:", paste(intersect(s1, s2), collapse = ", "), "\n")

# 与 D3 池交集
pool <- fread("D:/A生信分析/虚拟敲除/第3天执行结果/data/candidate_pool.csv")
p2m <- function(x) paste0(substr(x, 1, 1), tolower(substr(x, 2, nchar(x))))
pool_m <- unique(p2m(pool$gene))
for (g in c("Serpine1", "Tmsb4x")) {
  t100 <- head(ko2[target == g][order(p.value)]$gene, 100)
  cat(sprintf("%s: top100 命中 D3 池 %d 个: %s\n", g, sum(t100 %in% pool_m),
              paste(intersect(t100, pool_m), collapse = ", ")))
}
cat("\nD7-02d done\n")