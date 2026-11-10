# ============================================================
# D7-02b 虚拟敲除 v2：真正的网络扰动扩散
# v1 的问题：只删边 -> 除目标自身外所有基因 delta_degree 都是 1，无区分力
# v2 做法（贴近 scTenifoldKnk 思路）：
#   1. 用 PC regression 得到基因调控网络 A（基因 x 基因）
#   2. 敲除目标：把 A 中目标基因的入边与出边置 0
#   3. 在网络 A 上计算：敲除后每个基因到"核心调控模块"的可达性变化
#      —— 用重启随机游走（RWR）从目标基因出发，比较 A 与 A_ko 的稳态分布
#   4. delta = RWR_A - RWR_Ako，即"敲除导致的通路影响分数"
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(igraph)
})
out <- "D:/A生信分析/虚拟敲除/第7天执行结果"
A <- readRDS(file.path(out, "data", "D7_network_sym.rds"))
cat("network:", nrow(A), "genes\n")

# 保留非零边，转稀疏
A <- as.matrix(A)
A[abs(A) < 1e-8] <- 0
cat("nonzero edges:", sum(A != 0), "\n")

# ---- RWR（重启随机游走），返回稳态分布 ----
rwr <- function(W, seed_genes, restart = 0.15, tol = 1e-8, maxit = 500) {
  W <- abs(W)
  # 列归一化（列随机游走）
  cs <- colSums(W); cs[cs == 0] <- 1
  P <- sweep(W, 2, cs, "/")
  n <- nrow(P)
  p0 <- numeric(n); p0[seed_genes] <- 1 / length(seed_genes)
  p <- p0
  for (i in seq_len(maxit)) {
    p_new <- (1 - restart) * as.numeric(P %*% p) + restart * p0
    if (sum(abs(p_new - p)) < tol) { p <- p_new; break }
    p <- p_new
  }
  as.numeric(p)
}

run_ko_v2 <- function(target, topN = 30) {
  if (!target %in% rownames(A)) return(NULL)
  B <- A
  B[target, ] <- 0; B[, target] <- 0
  seeds <- which(colSums(A[, , drop = FALSE]) > 0)
  w0 <- rwr(A, seeds)
  w1 <- rwr(B, seeds)
  d  <- w0 - w1
  res <- data.table(gene = rownames(A), target = target,
                    rwr_orig = w0, rwr_ko = w1, delta = d)
  setorder(res, -delta)
  res
}

k1 <- run_ko_v2("Serpine1")
k2 <- run_ko_v2("Tmsb4x")
ko_all <- rbind(k1, k2)
fwrite(ko_all, file.path(out, "results", "D7_knockout_rwr.tsv"), sep = "\t")

cat("\n===== Serpine1 KO: top 20 受影响基因 =====\n")
print(k1[gene != "Serpine1"][1:20, .(gene, delta = signif(delta, 3))])
cat("\n===== Tmsb4x KO: top 20 受影响基因 =====\n")
print(k2[gene != "Tmsb4x"][1:20, .(gene, delta = signif(delta, 3))])

# ---- 与 D3 候选池交集 ----
pool <- fread("D:/A生信分析/虚拟敲除/第3天执行结果/data/candidate_pool.csv")
pool_m <- unique(paste0(substr(pool$gene, 1, 1), tolower(substr(pool$gene, 2, nchar(pool$gene)))))
for (nm in c("Serpine1", "Tmsb4x")) {
  kk <- ko_all[target == nm]
  top100 <- head(kk[gene != nm], 100)
  cat(sprintf("\n%s: top100 中落在 D3 候选池的 %d 个 (%s)\n", nm,
              sum(top100$gene %in% pool_m),
              paste(intersect(top100$gene, pool_m), collapse = ", ")))
}

# ---- 两个敲除的重叠 ----
ov <- merge(k1[, .(gene, d1 = delta)], k2[, .(gene, d2 = delta)], by = "gene")
ov[, agree := sign(d1) == sign(d2) & abs(d1) > 0 & abs(d2) > 0]
cat(sprintf("\n两敲除共同受影响基因 %d 个，方向一致 %d 个 (%.0f%%)\n",
            nrow(ov), sum(ov$agree), 100 * mean(ov$agree)))
fwrite(ov, file.path(out, "results", "D7_knockout_overlap.tsv"), sep = "\t")

# ---- 图 ----
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
mk <- function(dt, nm, n = 15) {
  d <- dt[gene != nm][order(-abs(delta))][1:n]
  d[, gene := factor(gene, levels = rev(gene))]
  ggplot(d, aes(x = gene, y = delta, fill = delta > 0)) +
    geom_col(width = 0.72) + coord_flip() +
    scale_fill_manual(values = c(`TRUE` = "#C44E52", `FALSE` = "#4C72B0"), guide = "none") +
    geom_hline(yintercept = 0, linewidth = 0.35) +
    labs(x = NULL, y = "RWR influence change (delta)",
         title = sprintf("%s virtual knockout", nm)) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0))
}
p <- mk(k1, "Serpine1") | mk(k2, "Tmsb4x")
ggsave(file.path(out, "figures", "D7_Fig2_knockout.png"), p, width = 13, height = 6.5, dpi = 300)
ggsave(file.path(out, "figures", "D7_Fig2_knockout.pdf"), p, width = 13, height = 6.5)

cat("\nD7-02b done\n")