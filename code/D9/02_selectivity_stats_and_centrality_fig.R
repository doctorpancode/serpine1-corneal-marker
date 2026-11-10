# ============================================================
# D9-02  待办 10.1-(2) + 10.1-(3)：Fig12 统计检验 + 中心性对照正式图
#        并额外给 10.2-(1) 的"真实敲低验证"留出可复用接口（D9-05 用）
#
#   10.1-(2)：D8 Fig12 只给了中位数与四分位，缺 p 值。
#             本脚本补：
#               - Wilcoxon 秩和检验（背景 vs 核心候选的 log2 选择性比值）
#               - 置换检验（10000 次，只重抽 11 个"核心"标签）
#               - Cliff's delta（非参数效应量）
#               - bootstrap 95% CI（中位数差）
#   10.1-(3)：中心性对照只存在于 tsv，图上没有。
#             本脚本画：
#               x = out_degree, y = 上皮 mean|Δ|
#               + Spearman rho=0.9169
#               + 11 个核心基因 ggrepel 标注
#               + 分层参考线（按出度三分位）
#
#   产出：
#     results/D9_fig12_stats.tsv
#     results/D9_centrality_by_stratum.tsv
#     figures/main/D9_Fig3_selectivity_test.{png,pdf}
#     figures/main/D9_Fig4_centrality_control.{png,pdf}
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggrepel); library(patchwork)
})
OUT <- "D:/A生信分析/虚拟敲除/第9天执行结果"
RES <- file.path(OUT, "results"); FIG <- file.path(OUT, "figures", "main")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

set.seed(20261110)

## ---------- 1. 读 Fig12 原始 data ----------
# 注意：D8 的结果都在"第8天执行结果/results"，不是第9天！
D8 <- "D:/A生信分析/虚拟敲除/第8天执行结果/results"
sel  <- fread(file.path(D8, "D8_epi_selectivity_ratio.tsv"))
imp  <- fread(file.path(D8, "D8_ko_impact_rank.tsv"))
edge <- fread(file.path(D8, "D8_network_edge_summary.tsv"))

CORE <- c("Tmsb4x","Serpine1","P4hb","S100a13","Sod1","Prdx2","Vdac1","Capg","Gja1","Xiap","Vim")
sel[, grp := factor(grp, levels = c("background", "core candidate"))]
cat("总 KO 数:", nrow(sel), " 核心:", sel[grp == "core candidate", .N],
    " 背景:", sel[grp == "background", .N], "\n")

x_bg <- sel[grp == "background", ratio]
x_cr <- sel[grp == "core candidate", ratio]

## ---------- 2. 检验：Wilcoxon + 置换 + Cliff's delta + bootstrap ----------
wl <- wilcox.test(x_cr, x_bg, exact = FALSE, alternative = "greater")
wl2 <- wilcox.test(x_cr, x_bg, exact = FALSE, alternative = "two.sided")

# 置换检验：固定 11 个"核心"标签随机分配到 947 个 KO 上
n_cr <- length(x_cr); n_all <- nrow(sel)
all_ratio <- sel$ratio
perm <- replicate(10000, {
  idx <- sample.int(n_all, n_cr)
  median(all_ratio[idx]) - median(all_ratio[-idx])
})
obs_diff <- median(x_cr) - median(x_bg)
perm_p <- (sum(perm >= obs_diff) + 1) / (length(perm) + 1)

# Cliff's delta：P(X>Y) - P(X<Y)
cd <- (sum(outer(x_cr, x_bg, ">")) - sum(outer(x_cr, x_bg, "<"))) / (length(x_cr) * length(x_bg))

# bootstrap 中位数差的 95% CI
boot <- replicate(5000, {
  median(sample(x_cr, replace = TRUE)) - median(sample(x_bg, replace = TRUE))
})
ci <- quantile(boot, c(0.025, 0.975))

stats <- data.table(
  test = c("wilcoxon_greater", "wilcoxon_two_sided", "permutation_10000",
           "cliffs_delta", "bootstrap_median_diff", "median_background",
           "median_core", "n_background", "n_core"),
  value = c(wl$p.value, wl2$p.value, perm_p, cd, obs_diff,
            median(x_bg), median(x_cr), length(x_bg), length(x_cr)),
  ci_lo = c(NA, NA, NA, NA, ci[1], NA, NA, NA, NA),
  ci_hi = c(NA, NA, NA, NA, ci[2], NA, NA, NA, NA)
)
fwrite(stats, file.path(RES, "D9_fig12_stats.tsv"), sep = "\t")
cat("\n=== Fig12 统计 ==="); print(stats)
cat("\nCliff's delta 解释：|d|<0.147 negligible, <0.33 small, <0.474 medium, 否则 large\n")

## ---------- 3. 中心性对照 ----------
# 出度直接从 WT 网络（947x947）算，不依赖只有 12 行的 edge_summary
WT <- readRDS("D:/A生信分析/虚拟敲除/第8天执行结果/data/D8_WT_network.rds")
M  <- as.matrix(WT); diag(M) <- 0
outdeg <- rowSums(M) != 0            # 逻辑：出边是否存在
outdeg_n <- rowSums(M != 0)          # 出度计数
indeg_n  <- colSums(M != 0)
deg_dt <- data.table(ko_gene = rownames(M),
                     out_degree = as.numeric(outdeg_n),
                     in_degree  = as.numeric(indeg_n))

net <- copy(imp)[mean_abs > 0]
net <- merge(net, deg_dt, by = "ko_gene")
# 上皮均 |Δ|（18 基因）
epi_eff <- sel[, .(ko_gene, epithel)]
net <- merge(net, epi_eff, by = "ko_gene")
for (cc in c("out_degree","in_degree","epithel","mean_abs")) storage.mode(net[[cc]]) <- "double"
cat("中心性对照可用 KO:", nrow(net), "\n")

rt <- suppressWarnings(cor.test(net$out_degree, net$epithel, method = "spearman", exact = FALSE))
cat("Spearman(out_degree, epithelial mean|delta|) rho =", round(unname(rt$estimate), 4),
    " p =", signif(rt$p.value, 3), "\n")

# 出度三分位分层后的组内相关（去掉中心性后还剩多少）
net[, strata := cut(out_degree, quantile(out_degree, c(0, 1/3, 2/3, 1)),
                    include.lowest = TRUE, labels = c("low", "mid", "high"))]
strata_tab <- net[, {
  r <- suppressWarnings(cor.test(out_degree, epithel, method = "spearman", exact = FALSE))
  .(n = .N, rho = unname(r$estimate), p = r$p.value,
    median_out = median(out_degree), median_epi = median(epithel))
}, by = strata]
fwrite(strata_tab, file.path(RES, "D9_centrality_by_stratum.tsv"), sep = "\t")
cat("\n=== 出度三分位分层 ===\n"); print(strata_tab)

# 残差化：用 out_degree 的 loess 拟合上皮效应，取残差排名
fit_res <- loess(epithel ~ out_degree, data = net, span = 0.5)
net[, epi_resid := resid(fit_res)]
net[, resid_rank := frank(-epi_resid, ties.method = "first")]
core_net <- net[ko_gene %in% CORE]
core_net[, core_rank_resid := match(ko_gene, net[order(-epi_resid), ko_gene])]
fwrite(core_net[order(-epi_resid),
                .(ko_gene, out_degree, in_degree, epithel, epi_resid, core_rank_resid)],
       file.path(RES, "D9_core_genes_residual_rank.tsv"), sep = "\t")
cat("\n=== 核心基因在去中心性残差下的排名 ===\n")
print(core_net[order(-epi_resid), .(ko_gene, out_degree, epithel = round(epithel, 5),
                                    epi_resid = round(epi_resid, 5), core_rank_resid)])

## ---------- 4. Fig：D9-Fig3 选择性检验 ----------
d <- copy(sel)
d[, y := as.numeric(grp)]
labs_bg <- unique(d[grp == "background", .(y = as.numeric(grp[1]), panel = "background")])
lab_dt <- d[grp == "core candidate", .(ko_gene, ratio,
                                       y = as.numeric(grp[1]))]
p3 <- ggplot(d, aes(x = ratio, y = y)) +
  geom_violin(aes(fill = grp), scale = "width", linewidth = 0.25, alpha = 0.55) +
  geom_boxplot(width = 0.22, outlier.shape = NA, linewidth = 0.3, fill = "white") +
  geom_jitter(data = d[grp == "background"], width = 0, height = 0.10,
              size = 0.35, alpha = 0.30, colour = "grey45") +
  geom_jitter(data = d[grp == "core candidate"], width = 0, height = 0.075,
              size = 1.5, alpha = 0.9, colour = "#B2182B") +
  geom_segment(data = lab_dt, aes(x = ratio, xend = ratio, y = y, yend = y + 0.42),
               linewidth = 0.2, colour = "#B2182B") +
  geom_text_repel(data = lab_dt, aes(x = ratio, y = y + 0.44, label = ko_gene),
                  inherit.aes = FALSE, size = 2.2,
                  min.segment.length = Inf, direction = "both",
                  box.padding = 0.22, point.padding = 0, max.overlaps = 30,
                  segment.colour = NA, colour = "#B2182B") +
  scale_fill_manual(values = c(background = "grey78", `core candidate` = "#F4A582")) +
  scale_y_continuous(breaks = 1:2, labels = c("background", "core candidate"),
                     limits = c(0.5, 2.9), expand = c(0, 0)) +
  labs(x = expression(log[2]~"(epithelial mean"~"|"*Delta*"| / other mean"~"|"*Delta*"|)"),
       y = NULL,
       title = "Epithelial selectivity of candidate genes vs the whole network",
       subtitle = sprintf("median %.2f (n = %d) vs %.2f (n = %d)\nWilcoxon p = %.2f  |  Cliff's delta = %.2f (negligible)  |  permutation p = %.2f",
                          median(x_bg), length(x_bg), median(x_cr), length(x_cr),
                          wl$p.value, cd, perm_p)) +
  theme_bw(base_size = 9) +
  theme(legend.position = "none", panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.subtitle = element_text(size = 6.6, colour = "grey30"))
ggsave(file.path(FIG, "D9_Fig3_selectivity_test.png"), p3, width = 6.6, height = 3.2, dpi = 300)
ggsave(file.path(FIG, "D9_Fig3_selectivity_test.pdf"), p3, width = 6.6, height = 3.2)
cat("\n[saved] D9_Fig3_selectivity_test\n")

## ---------- 5. Fig：D9-Fig4 中心性对照 ----------
p4 <- ggplot(net, aes(x = out_degree, y = epithel)) +
  geom_point(size = 0.55, alpha = 0.30, colour = "grey35") +
  geom_smooth(method = "loess", span = 0.6, se = TRUE, linewidth = 0.6,
              colour = "#2166AC", fill = "grey85") +
  geom_point(data = core_net, aes(x = out_degree, y = epithel),
             size = 2.0, shape = 21, fill = "#B2182B", colour = "black", stroke = 0.3) +
  geom_text_repel(data = core_net, aes(x = out_degree, y = epithel, label = ko_gene),
                  size = 2.4, box.padding = 0.35, point.padding = 0.15,
                  min.segment.length = 0, segment.size = 0.2, max.overlaps = 30) +
  labs(x = "Network out-degree (transcription factor-like potential)",
       y = expression("Epithelial effect  mean"~"|"*Delta*"| (18 genes)"),
       title = "Perturbation magnitude tracks network centrality",
       subtitle = sprintf("n = %d knockouts   Spearman rho = %.4f, p < 1e-300   core candidates highlighted",
                          nrow(net), unname(rt$estimate))) +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = 7, colour = "grey30"))
ggsave(file.path(FIG, "D9_Fig4_centrality_control.png"), p4, width = 5.4, height = 4.0, dpi = 300)
ggsave(file.path(FIG, "D9_Fig4_centrality_control.pdf"), p4, width = 5.4, height = 4.0)
cat("[saved] D9_Fig4_centrality_control\n")

## ---------- 6. 拼一张双面板（Fig6 直接用） ----------
pcombo <- p4 + p3 + plot_layout(widths = c(1.15, 1)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 10, face = "bold"))
ggsave(file.path(FIG, "D9_Fig5_selectivity_and_centrality_combo.png"), pcombo,
       width = 11.6, height = 4.0, dpi = 300)
ggsave(file.path(FIG, "D9_Fig5_selectivity_and_centrality_combo.pdf"), pcombo,
       width = 11.6, height = 4.0)
cat("[saved] D9_Fig5_selectivity_and_centrality_combo\n")

cat("\n=== DONE 02 ===\n")