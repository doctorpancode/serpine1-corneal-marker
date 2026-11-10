# ============================================================
# D7-04b 主图更新：Fig3 scorecard / Fig8 虚拟敲除（官方 scTenifoldKnk 版）
# 只重画受 02d 影响的图，其它 7 张沿用 04_main_figures.R 产物
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})
D7 <- "D:/A生信分析/虚拟敲除/第7天执行结果"
D6 <- "D:/A生信分析/虚拟敲除/第6天执行结果"
D5 <- "D:/A生信分析/虚拟敲除/第5天执行结果"
D4 <- "D:/A生信分析/虚拟敲除/第4天执行结果"
fig <- file.path(D7, "figures", "main")

ko  <- fread(file.path(D7, "results", "D7_knk_official_twoKO.tsv"))
gs  <- fread(file.path(D5, "results", "D5_group_stats.tsv"))
sc  <- fread(file.path(D4, "results", "D4_candidate_scorecard.tsv"))
cat("KO rows:", nrow(ko), "\n")

# ---------- Fig8：两个敲除的 Z 值条形 ----------
mk8 <- function(nm, n = 12) {
  d <- ko[target == nm & gene != nm][order(-Z)][1:n]
  d[, gene := factor(gene, levels = rev(gene))]
  d[, sig := p.adj < 0.05]
  ggplot(d, aes(x = gene, y = Z, fill = sig)) +
    geom_col(width = 0.72) + coord_flip() +
    geom_hline(yintercept = 1.96, linetype = 2, linewidth = 0.4, colour = "grey40") +
    scale_fill_manual(values = c(`TRUE` = "#C44E52", `FALSE` = "#B8C4CE"),
                      labels = c("TRUE" = "FDR < 0.05", "FALSE" = "FDR ≥ 0.05"),
                      name = NULL) +
    labs(x = NULL, y = "differential-regulation Z",
         title = sprintf("%s virtual knockout (scTenifoldKnk)", nm)) +
    theme_bw(base_size = 9.5) +
    theme(plot.title = element_text(face = "bold", size = 9.5),
          legend.position = "bottom", legend.key.size = unit(0.3, "cm"))
}
p8 <- mk8("Serpine1") | mk8("Tmsb4x")
ggsave(file.path(fig, "Fig8_virtual_knockout.png"), p8, width = 11, height = 4.8, dpi = 300)
ggsave(file.path(fig, "Fig8_virtual_knockout.pdf"), p8, width = 11, height = 4.8)

# ---------- Fig3：证据记分卡（加官方敲除命中列）----------
pool <- fread("D:/A生信分析/虚拟敲除/第3天执行结果/data/candidate_pool.csv")
p2m <- function(x) paste0(substr(x, 1, 1), tolower(substr(x, 2, nchar(x))))
pool[, mgene := p2m(gene)]
top50 <- list()
for (g in c("Serpine1", "Tmsb4x")) top50[[g]] <- ko[target == g & gene != g][order(p.value)][1:50]$gene

s <- copy(sc)
cn <- setdiff(names(s), c("gene"))
if (!"KO_hit_in_top50" %in% cn) {
  s[, KO_hit_in_top50 := sapply(gene, function(x) {
    xm <- p2m(x)
    hits <- c(if (xm %in% top50$Serpine1) "Serpine1", if (xm %in% top50$Tmsb4x) "Tmsb4x")
    if (length(hits)) paste(hits, collapse = "+") else "none"
  })]
}
cat("scorecard cols:", paste(names(s), collapse = ", "), "\n")
print(s[, .N, by = KO_hit_in_top50])

# 第4天 scorecard 没有 ml_rank 列，用 score 排序（score 越大证据越强）
s[, rk := frank(-score, ties.method = "min")]
lab <- s[order(rk)][1:min(15, .N)]
lab[, gene := factor(gene, levels = rev(gene))]
lab[, ko_lab := fifelse(KO_hit_in_top50 != "none", "hit", "no hit")]
lab[, ko_lab := factor(ko_lab, levels = c("hit", "no hit"))]
p3 <- ggplot(lab, aes(y = gene, x = score, fill = ko_lab)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = score), hjust = -0.25, size = 3, colour = "grey20") +
  scale_fill_manual(values = c(hit = "#C44E52", `no hit` = "#B8C4CE"),
                    name = "virtual KO\n(top-50)", drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(x = "evidence score (D4 scorecard)", y = NULL,
       title = "Candidate evidence scorecard with virtual-knockout annotation") +
  theme_bw(base_size = 9.5) +
  theme(plot.title = element_text(face = "bold", size = 9.5),
        legend.key.size = unit(0.3, "cm"))
ggsave(file.path(fig, "Fig3_scorecard.png"), p3, width = 6.6, height = 4.6, dpi = 300)
ggsave(file.path(fig, "Fig3_scorecard.pdf"), p3, width = 6.6, height = 4.6)

# ---------- Fig1 文案：把 scTenifoldKnk 改成官方 ----------
cat("\nFig8/Fig3 已按官方 scTenifoldKnk 结果重画\n")
cat("Serpine1 KO top50:", paste(top50$Serpine1, collapse = ", "), "\n")
cat("Tmsb4x KO top50:", paste(top50$Tmsb4x, collapse = ", "), "\n")
cat("\nD7-04b done\n")