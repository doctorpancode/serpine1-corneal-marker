# ============================================================
# D8-04b 修版：Fig10 影响谱热图重绘 + Fig12 改为散点对比
#   问题1：Fig10 原版 z-score 在三值尺度上把 947x120 压成白纸，
#          只有 core 行有信号 → 改用行内归一（每个 KO 谱按自身最大值缩放）
#          + 只保留信号最强的列子集，并对绝对值取 log10 着色。
#   问题2：Fig12 原版 947 个 x 轴刻度标签全糊 → 改为
#          x = mean|delta|(全部基因), y = mean|delta|(上皮程序)，散点 + 对角参考线。
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(ComplexHeatmap); library(circlize); library(grid); library(svglite)
})
out <- "D:/A生信分析/虚拟敲除/第8天执行结果"
fig <- file.path(out, "figures", "main")

spec <- fread(file.path(out, "results", "D8_knockdown_spectrum.tsv"))
imp  <- fread(file.path(out, "results", "D8_ko_impact_rank.tsv"))
raw  <- readRDS(file.path(out, "data", "D8_scan_raw.rds"))

CORE <- c("Serpine1", "Tmsb4x", "P4hb", "S100a13", "Sod1", "Prdx2",
          "Vdac1", "Capg", "Gja1", "Xiap", "Vim")
CORE <- intersect(CORE, unique(spec$ko_gene))
EPI  <- c("Krt12", "Aqp5", "Gsta4", "Upk3bl", "Upk1b", "Krt15", "Krt17",
          "Krt19", "Krt13", "Dsg1a", "Perp", "Tmprss11d", "Lypd3", "Calml3",
          "Gsta2", "Sfn", "Krt5", "Krt14")
EPI  <- intersect(EPI, unique(spec$affected_gene))
cat("core KOs:", length(CORE), "| epi genes:", length(EPI), "\n")

save_all <- function(p, name, w, h) {
  ggsave(file.path(fig, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  ggsave(file.path(fig, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(fig, paste0(name, ".svg")), p, width = w, height = h, device = svglite::svglite)
}

# ---------- Fig10 重绘 ----------
D <- abs(raw$perturbationDistances)
D[is.na(D)] <- 0
# 每个 KO 谱按自身第 95 百分位缩放 → 移除"敲除强度"造成的行间量级差，
# 只留"打在哪里"的形状差异。
scale_row <- function(v) {
  q <- quantile(v[v > 0], 0.95, na.rm = TRUE)
  if (!is.finite(q) || q <= 0) q <- max(v, na.rm = TRUE)
  if (!is.finite(q) || q <= 0) return(rep(0, length(v)))
  v / q
}
Dn <- t(apply(D, 1, scale_row))
dimnames(Dn) <- dimnames(D)

# 选列：上皮程序 + 每行 top5 命中频次最高的基因，取前 90
freq <- sort(table(spec[order(ko_gene, -abs_distance)][, head(.SD, 5), by = ko_gene]$affected_gene),
              decreasing = TRUE)
# 只保留在 >=10% 敲除里都有信号的基因，否则会带进一整片空白列
prev <- colMeans(Dn > 0.05)
col_genes <- unique(c(EPI, names(freq)))
col_genes <- intersect(col_genes, colnames(Dn))
col_genes <- col_genes[prev[col_genes] >= 0.10]
col_genes <- head(col_genes, 80)
cat("kept columns:", length(col_genes), "| epithelial in panel:", sum(col_genes %in% EPI), "\n")

row_ord <- c(CORE, setdiff(imp[order(-mean_abs)]$ko_gene, CORE))
row_ord <- row_ord[row_ord %in% rownames(Dn)]
M <- Dn[row_ord, col_genes, drop = FALSE]
M[M > 1.5] <- 1.5   # 顶部截断，避免少数极值吃掉动态范围

col_ha <- HeatmapAnnotation(
  program = ifelse(colnames(M) %in% EPI, "epithelial", "other"),
  col = list(program = c(epithelial = "#C44E52", other = "#B8C4CE")),
  show_annotation_name = FALSE,
  simple_anno_size = unit(3, "mm"),
  annotation_legend_param = list(title = "affected gene", title_gp = gpar(fontsize = 8),
                                 labels_gp = gpar(fontsize = 7))
)
row_ha <- rowAnnotation(
  core = anno_mark(at = match(CORE, rownames(M)), labels = CORE,
                   labels_gp = gpar(fontsize = 7.5, fontface = "bold"),
                   link_width = unit(7, "mm"), padding = unit(1.4, "mm")),
  width = unit(24, "mm")
)
ht <- Heatmap(M, name = "scaled |delta|\n(row-normalised)",
              col = colorRamp2(c(0, 0.25, 0.6, 1.0, 1.5),
                               c("#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#8B1A1A")),
              top_annotation = col_ha, right_annotation = row_ha,
              cluster_rows = FALSE, cluster_columns = TRUE,
              show_row_names = FALSE, show_column_names = TRUE,
              column_names_gp = gpar(fontsize = 5.2),
              column_names_rot = 90,
              row_split = factor(ifelse(rownames(M) %in% CORE, "core candidates", "background"),
                                 levels = c("core candidates", "background")),
              row_title_gp = gpar(fontsize = 8),
              column_title = sprintf("Knockout impact spectrum: %d virtual knockouts x %d readout genes",
                                     nrow(M), ncol(M)),
              column_title_gp = gpar(fontsize = 10, fontface = "bold"),
              heatmap_legend_param = list(title_gp = gpar(fontsize = 7.5),
                                          labels_gp = gpar(fontsize = 7),
                                          legend_height = unit(30, "mm")))

for (dev in c("pdf", "png")) {
  if (dev == "pdf") pdf(file.path(fig, "Fig10_knockout_impact_spectrum.pdf"), width = 10, height = 6.6)
  else png(file.path(fig, "Fig10_knockout_impact_spectrum.png"), width = 10, height = 6.6,
           units = "in", res = 300)
  draw(ht, merge_legend = TRUE); dev.off()
}
svglite::svglite(file.path(fig, "Fig10_knockout_impact_spectrum.svg"), width = 10, height = 6.6)
draw(ht, merge_legend = TRUE); dev.off()
cat("Fig10 redrawn\n")

# ---------- Fig12 重做：选择性比值分布 ----------
cmp <- spec[!is.na(abs_distance), .(
  epithel = mean(abs_distance[affected_gene %in% EPI]),
  other   = mean(abs_distance[!affected_gene %in% EPI])
), by = ko_gene]
cmp <- merge(cmp, imp[, .(ko_gene, mean_abs)], by = "ko_gene")
cmp[, core := ko_gene %in% CORE]
cmp <- cmp[epithel > 0 & other > 0]
cmp[, ratio := log2(epithel / other)]
cmp[, grp := ifelse(core, "core candidate", "background")]
cmp[, grp := factor(grp, levels = c("background", "core candidate"))]
cat("ratio (log2) summary by group:\n"); print(cmp[, .(n = .N, median = median(ratio),
                                                     q25 = quantile(ratio, .25),
                                                     q75 = quantile(ratio, .75)), by = grp])
cat("\ncore rows by ratio:\n"); print(cmp[core == TRUE][order(-ratio), .(ko_gene, epithel, other, ratio)])

p12 <- ggplot(cmp, aes(x = ratio, y = grp)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_jitter(data = cmp[grp == "background"], height = 0.08, colour = "grey70",
              size = 0.9, alpha = 0.55) +
  geom_boxplot(data = cmp[grp == "background"], width = 0.45, outlier.shape = NA,
               fill = NA, colour = "grey35", linewidth = 0.45) +
  geom_point(data = cmp[grp == "core candidate"], aes(colour = ko_gene), size = 3.0) +
  ggrepel::geom_text_repel(data = cmp[grp == "core candidate"], aes(label = ko_gene, colour = ko_gene),
                           size = 3, fontface = "bold", show.legend = FALSE,
                           nudge_y = 0.42, direction = "both", box.padding = 0.55, force = 3,
                           min.segment.length = 0.1, max.overlaps = 20, segment.size = 0.3) +
  scale_colour_manual(values = setNames(
    colorRampPalette(c("#B2182B", "#762A83", "#2166AC"))(nrow(cmp[core == TRUE])),
    cmp[core == TRUE][order(-ratio)]$ko_gene), guide = "none") +
  labs(x = sprintf("log2 [ effect on %d-gene epithelial program / effect on all other genes ]", length(EPI)),
       y = NULL,
       title = "Epithelial selectivity of virtual knockouts",
       subtitle = sprintf("positive = the knockout hits epithelial genes harder than the rest of the network. %d knockouts; box = all non-core genes",
                          nrow(cmp))) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 10.5),
        plot.subtitle = element_text(size = 7.6, colour = "grey40"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 9))
save_all(p12, "Fig12_epithelial_vs_rest", 8.4, 2.9)
cat("Fig12 redrawn\n")

fwrite(cmp[order(-ratio)], file.path(out, "results", "D8_epi_selectivity_ratio.tsv"), sep = "\t")
cat("\nD8-04b done\n")

