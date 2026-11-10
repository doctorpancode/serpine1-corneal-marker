# ============================================================
# D8-03 第三个敲除：P4HB（以及 S100A13 / VIM 作为并列验证）
#
# 目的：D7 只敲了 Serpine1 / Tmsb4x，两者都扰动"上皮角化程序"。
#       需要第三个独立基因来回答："上皮程序被扰动"是这两个基因的
#       特性，还是任何核心基因敲了都会看到（即网络中心性伪影）。
#
# 方法：与 D7-02d 完全同构——同一张 WT 网络、同一份 qcCounts、
#       官方 manifoldAlignment + scTenifoldKnk::dRegulation。
#       区别只在网络规模从 494 涨到 D8 的 900+，因此同时给出
#       "同网络内部可比"（Serpine1/Tmsb4x/P4hb 一起跑）与
#       "跨网络口径"（与 D7 的 494 网对比）两套读数。
#
# 输出：
#   results/D8_third_ko.tsv             多敲除合并长表
#   results/D8_third_ko_summary.tsv     每个敲除的显著基因数 / top 基因 / 方向
#   results/D8_sig_overlap.tsv          显著基因两两重叠矩阵（含 D7 口径）
#   results/D8_epithelial_program_score.tsv  上皮程序被扰动程度打分
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(scTenifoldKnk); library(scTenifoldNet)
})
set.seed(20261110)
out <- "D:/A生信分析/虚拟敲除/第8天执行结果"

WT  <- readRDS(file.path(out, "data", "D8_WT_network.rds"))
qcX <- readRDS(file.path(out, "data", "D8_knk_qcCounts.rds"))
stopifnot(identical(rownames(WT), rownames(qcX)))

KO_LIST <- c("Serpine1", "Tmsb4x", "P4hb", "S100a13", "Vim")
KO_LIST <- intersect(KO_LIST, rownames(WT))
cat("knockouts to run:", paste(KO_LIST, collapse = ", "), "\n")

# ---------------- 1. 逐个官方流程 ----------------
allres <- list()
for (g in KO_LIST) {
  cat("\n===== knocking out", g, "=====\n")
  KO <- WT
  KO[g, ] <- 0
  cat("  WT outdegree:", sum(WT[g, ] != 0), "-> KO outdegree:", sum(KO[g, ] != 0), "\n")

  set.seed(20261110)
  dir <- tryCatch(
    scTenifoldKnk::knockoutDirection(qcX, gKO = list(g), genes = rownames(WT), t = 5),
    error = function(e) { cat("  knockoutDirection failed:", conditionMessage(e), "\n"); NULL })
  if (!is.null(dir)) dir <- setNames(as.numeric(dir), colnames(dir))

  set.seed(20261110)
  MA <- scTenifoldNet::manifoldAlignment(WT, KO, d = 2, nCores = 6)
  DR <- scTenifoldKnk::dRegulation(manifoldOutput = MA, gKO = g,
                                   empiricalNull = TRUE, direction = dir)
  dr <- as.data.frame(DR, stringsAsFactors = FALSE)
  if (!"gene" %in% names(dr)) dr[["gene"]] <- rownames(DR)
  if (!"p.adj" %in% names(dr)) dr[["p.adj"]] <- dr[["pvalue"]]
  if (!"p.value" %in% names(dr)) dr[["p.value"]] <- dr[["p.adj"]]
  dr[["target"]] <- g
  dr[["is_target"]] <- dr[["gene"]] == g
  cat("  列名:", paste(names(dr), collapse = ","), "\n")
  cat("  significant (FDR<0.05):", sum(dr$p.adj < 0.05, na.rm = TRUE),
      "| non-target:", sum(dr$p.adj < 0.05 & dr$gene != g, na.rm = TRUE), "\n")
  dr <- as.data.table(dr)
  nn <- min(10, nrow(dr[gene != g]))
  if (nn > 0) print(dr[gene != g][order(p.value)][1:nn])
  allres[[g]] <- dr
}
ko <- rbindlist(allres, fill = TRUE)
fwrite(ko, file.path(out, "results", "D8_third_ko.tsv"), sep = "\t")

# ---------------- 2. 摘要 ----------------
summ <- ko[, .(
  n_sig_total     = sum(p.adj < 0.05, na.rm = TRUE),
  n_sig_nontarget = sum(p.adj < 0.05 & !is_target, na.rm = TRUE),
  top5_nontarget  = paste(head(gene[!is_target][order(p.value[!is_target])], 5), collapse = ","),
  n_up            = sum(p.adj < 0.05 & !is_target & direction == "up", na.rm = TRUE),
  n_down          = sum(p.adj < 0.05 & !is_target & direction == "down", na.rm = TRUE),
  outdegree       = sum(WT[gene[1], ] != 0),
  indegree        = sum(WT[, gene[1]] != 0)
), by = target]
setorder(summ, -n_sig_nontarget)
print(summ)
fwrite(summ, file.path(out, "results", "D8_third_ko_summary.tsv"), sep = "\t")

# ---------------- 3. 显著基因两两重叠 ----------------
sig <- lapply(KO_LIST, function(g) ko[target == g & p.adj < 0.05 & !is_target & !is.na(p.adj)]$gene)
names(sig) <- KO_LIST
ov <- as.data.table(expand.grid(a = KO_LIST, b = KO_LIST, stringsAsFactors = FALSE))
ov[, n_overlap := mapply(function(x, y) length(intersect(sig[[x]], sig[[y]])), a, b)]
ov[, jaccard := mapply(function(x, y) {
  u <- union(sig[[x]], sig[[y]]); if (length(u) == 0) NA_real_ else length(intersect(sig[[x]], sig[[y]]))/length(u)
}, a, b)]
ov[, list_a := sapply(a, function(x) length(sig[[x]]))]
print(dcast(ov, a ~ b, value.var = "n_overlap"))
fwrite(ov, file.path(out, "results", "D8_sig_overlap.tsv"), sep = "\t")
cat("\n共享的显著扰动基因：\n")
tab <- table(unlist(sig))
print(sort(tab, decreasing = TRUE))

# ---------------- 4. 上皮程序打分 ----------------
# 上皮/角化程序基因集：D7 两敲共同显著的 3 个 + 富集 top 里的角蛋白/分化基因
EPI <- c("Krt12", "Aqp5", "Gsta4", "Upk3bl", "Upk1b", "Krt15", "Krt17",
         "Krt19", "Krt13", "Dsg1a", "Perp", "Tmprss11d", "Lypd3", "Calml3",
         "Gsta2", "Sfn", "Krt5", "Krt14")
EPI <- intersect(EPI, rownames(WT))
cat("\nepithelial gene set size in network:", length(EPI), "\n")

epi_score <- ko[gene %in% EPI & !is_target, .(
  n_epi_tested     = .N,
  n_epi_sig        = sum(p.adj < 0.05, na.rm = TRUE),
  mean_Z_epi       = mean(Z, na.rm = TRUE),
  median_p_epi     = median(p.value, na.rm = TRUE),
  n_epi_up         = sum(p.adj < 0.05 & direction == "up", na.rm = TRUE),
  n_epi_down       = sum(p.adj < 0.05 & direction == "down", na.rm = TRUE),
  top_epi_hits     = paste(head(gene[order(p.value)], 5), collapse = ",")
), by = target]
setorder(epi_score, -mean_Z_epi)
print(epi_score)
fwrite(epi_score, file.path(out, "results", "D8_epithelial_program_score.tsv"), sep = "\t")

# ---------------- 5. 网络拓扑阴性对照 ----------------
# 如果"敲谁都扰上皮"，那上皮扰动强度应与目标基因的拓扑中心性相关。
# 用 02 的全基因扫描结果做背景，算出上皮程序被扰动强度与
# "该基因出度"的相关性，作为中心性伪影的量化证据。
if (file.exists(file.path(out, "results", "D8_knockdown_spectrum.tsv"))) {
  spec <- fread(file.path(out, "results", "D8_knockdown_spectrum.tsv"))
  bg <- spec[affected_gene %in% EPI, .(mean_abs_epi = mean(abs_distance, na.rm = TRUE)), by = ko_gene]
  bg[, out_deg := rowSums(WT[ko_gene, , drop = FALSE] != 0)]
  ct <- suppressWarnings(cor.test(bg$mean_abs_epi, bg$out_deg, method = "spearman"))
  cat(sprintf("\n[中心性对照] 全 %d 个敲除 x 上皮基因集：rho(mean|dist|, outdegree) = %.3f, p = %.2e\n",
              nrow(bg), ct$estimate, ct$p.value))
  fwrite(data.table(metric = c("n_ko", "spearman_rho", "p_value"),
                    value = c(nrow(bg), round(unname(ct$estimate), 4), signif(ct$p.value, 3))),
         file.path(out, "results", "D8_centrality_null_check.tsv"), sep = "\t")
  bg[, is_core := ko_gene %in% c("Serpine1", "Tmsb4x", "P4hb", "S100a13", "Vim")]
  print(bg[is_core == TRUE][order(-mean_abs_epi)])
}

cat("\nD8-03 done\n")