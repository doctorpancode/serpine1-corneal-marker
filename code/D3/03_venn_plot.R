# D3-03 三主题 Venn 图 + 候选池分层图
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
out <- "D:/A生信分析/虚拟敲除/第3天执行结果"

venn <- fread(file.path(out, "results", "D3_venn_counts.tsv"))
pool <- fread(file.path(out, "data", "candidate_pool_full.csv"))
print(venn)

# ---- 手工精确绘制三圆 Venn（面积按真实计数缩放） ----
A <- venn[region == "A", n_genes]; B <- venn[region == "B", n_genes]; C <- venn[region == "C", n_genes]
AB <- venn[region == "A&B", n_genes]; AC <- venn[region == "A&C", n_genes]
BC <- venn[region == "B&C", n_genes]; ABC <- venn[region == "A&B&C", n_genes]
a_only <- A - AB - AC + ABC
b_only <- B - AB - BC + ABC
c_only <- C - AC - BC + ABC
ab_only <- AB - ABC; ac_only <- AC - ABC; bc_only <- BC - ABC

circ <- function(cx, cy, r, fill, alpha = 0.35) {
  t <- seq(0, 2 * pi, length.out = 200)
  data.frame(x = cx + r * cos(t), y = cy + r * sin(t), g = paste0(cx, cy))
}
d1 <- circ(-0.55, 0.35, 1, "#4C72B0")
d2 <- circ(0.55, 0.35, 1, "#55A868")
d3 <- circ(0, -0.55, 1, "#C44E52")
labs <- data.table(
  lab = c(a_only, b_only, c_only, ab_only, ac_only, bc_only, ABC),
  x = c(-1.05, 1.05, 0.0, 0.0, -0.72, 0.72, 0.0),
  y = c(0.55, 0.55, -1.0, 0.62, -0.35, -0.35, 0.12)
)
p1 <- ggplot() +
  geom_polygon(data = d1, aes(x, y), fill = "#4C72B0", alpha = 0.28) +
  geom_polygon(data = d2, aes(x, y), fill = "#55A868", alpha = 0.28) +
  geom_polygon(data = d3, aes(x, y), fill = "#C44E52", alpha = 0.28) +
  geom_text(data = labs, aes(x, y, label = lab), size = 4.2, fontface = "bold") +
  annotate("text", x = -1.5, y = 1.55, label = "A  CNV direct\n(CTD+DisGeNET+DISEASES2025+HPO)",
           size = 3.2, hjust = 0) +
  annotate("text", x = 0.45, y = 1.55, label = "B  Pathological neovascularization\n∩ corneal context",
           size = 3.2, hjust = 0) +
  annotate("text", x = -0.55, y = -2.0, label = "C  Corneal / limbal injury\n∩ WGCNA module genes",
           size = 3.2, hjust = 0) +
  coord_equal(xlim = c(-2.1, 2.4), ylim = c(-2.2, 1.9)) +
  theme_void(base_size = 12) +
  ggtitle("Three-theme evidence overlap (D3)")

# ---- 分层条形图 ----
pool[, tier := fifelse(in_theme_A + in_theme_B + in_theme_C == 3 & (in_GSE241715 == 1),
                       "triple ∧ detected", "other")]
tt <- fread(file.path(out, "results", "D3_pool_tiering.tsv"))
tt[, tier := factor(tier, levels = c("T1_triple_DEG", "T2_multi_DEG", "T3_DEG",
                                     "T4_triple_notDEG", "T5_other"))]
p2 <- ggplot(tt, aes(tier, N, fill = tier)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = N), vjust = -0.4, size = 4) +
  scale_fill_manual(values = c("#C44E52", "#DD8452", "#4C72B0", "#8172B3", "grey75")) +
  labs(x = NULL, y = "Number of genes",
       title = "Candidate pool tiering",
       subtitle = "T1 triple-theme AND DEG = 36 ; final pool (all DEG) = 167") +
  theme_bw(base_size = 12) + theme(legend.position = "none",
                                   axis.text.x = element_text(angle = 18, hjust = 1))

p <- p1 / p2 + plot_layout(heights = c(2.1, 1))
ggsave(file.path(out, "figures", "D3_Fig1_venn.pdf"), p, width = 7.6, height = 9.4)
ggsave(file.path(out, "figures", "D3_Fig1_venn.png"), p, width = 7.6, height = 9.4, dpi = 300)
cat("done\n")