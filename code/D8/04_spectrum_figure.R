# ============================================================
# D8-04 敲除影响谱图（Fig10 / Fig11 / Fig12） + 第三敲除图（Fig13）
#
# 输入：
#   results/D8_knockdown_spectrum.tsv          长表 ko x affected
#   results/D8_perturbation_distance_matrix.tsv 宽表
#   results/D8_ko_impact_rank.tsv              敲除强度排名
#   results/D8_third_ko.tsv                    第三敲除差异调控
#   results/D8_epithelial_program_score.tsv
#   results/D8_sig_overlap.tsv
#   data/D8_scan_raw.rds                       官方原始返回（给 perturbationMap）
#
# 输出（figures/main/，每张 PDF + PNG 300dpi + SVG 矢量）：
#   Fig10_knockout_impact_spectrum   全基因敲除影响谱热图（核心基因高亮）
#   Fig11_impact_landscape           PCA/扰动图式二维景观
#   Fig12_epithelial_vs_rest        上皮程序 vs 其余基因的效应对比
#   Fig13_third_knockout             P4HB/S100A13/VIM 第三敲除验证
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(ComplexHeatmap); library(circlize); library(grid); library(svglite)
})
out <- "D:/A生信分析/虚拟敲除/第8天执行结果"
fig <- file.path(out, "figures", "main")
dir.create(fig, showWarnings = FALSE, recursive = TRUE)

spec  <- fread(file.path(out, "results", "D8_knockdown_spectrum.tsv"))
imp   <- fread(file.path(out, "results", "D8_ko_impact_rank.tsv"))
ko3   <- fread(file.path(out, "results", "D8_third_ko.tsv"))
epi   <- fread(file.path(out, "results", "D8_epithelial_program_score.tsv"))
raw   <- readRDS(file.path(out, "data", "D8_scan_raw.rds"))

CORE <- c("Serpine1", "Tmsb4x", "P4hb", "S100a13", "Sod1", "Prdx2",
          "Vdac1", "Capg", "Gja1", "Xiap", "Vim")
CORE <- intersect(CORE, unique(spec$ko_gene))
EPI  <- c("Krt12", "Aqp5", "Gsta4", "Upk3bl", "Upk1b", "Krt15", "Krt17",
          "Krt19", "Krt13", "Dsg1a", "Perp", "Tmprss11d", "Lypd3", "Calml3",
          "Gsta2", "Sfn", "Krt5", "Krt14")
EPI  <- intersect(EPI, unique(spec$affected_gene))
cat("core KOs:", paste(CORE, collapse = ","), "\nepi genes:", length(EPI), "\n")

save_all <- function(p, name, w, h) {
  ggsave(file.path(fig, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  ggsave(file.path(fig, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(fig, paste0(name, ".svg")), p, width = w, height = h, device = svglite::svglite)
}

# ---------- Fig10 影响谱热图 ----------
# 取每个敲除谱里距离最大的前 40 个受影响基因的并集 + 上皮程序基因，
# 行 = 敲除（核心基因在最上），列 = 受影响基因，值 = |distance| 的列内 z 分数。
top_per_ko <- spec[order(ko_gene, -abs_distance)][, head(.SD, 40), by = ko_gene]
col_genes <- unique(c(top_per_ko$affected_gene, EPI))
col_genes <- head(col_genes, 120)

M <- raw$perturbationDistances[rownames(raw$perturbationDistances) %in% imp$ko_gene,
                               colnames(raw$perturbationDistances) %in% col_genes, drop = FALSE]
M <- abs(M)
M[is.na(M)] <- 0
# 列内 z-score（让"哪些敲除对这一列最狠"可比）
Mz <- scale(M)
Mz[!is.finite(Mz)] <- 0
Mz <- pmin(pmax(Mz, -3), 3)

row_ord <- c(CORE, setdiff(imp[order(-mean_abs)]$ko_gene, CORE))
row_ord <- row_ord[row_ord %in% rownames(Mz)]
Mz <- Mz[row_ord, , drop = FALSE]

col_ha <- HeatmapAnnotation(
  program = ifelse(colnames(Mz) %in% EPI, "epithelial/keratinization", "other"),
  col = list(program = c("epithelial/keratinization" = "#C44E52", "other" = "#B8C4CE")),
  show_annotation_name = FALSE,
  annotation_legend_param = list(title = "affected gene program")
)
row_ha <- rowAnnotation(
  core = anno_mark(at = match(CORE, rownames(Mz)), labels = CORE,
                   labels_gp = gpar(fontsize = 7, fontface = "bold"),
                   link_width = unit(6, "mm"), padding = unit(1.2, "mm")),
  width = unit(22, "mm")
)
ht <- Heatmap(Mz, name = "|delta| z-score\n(within gene)",
              col = colorRamp2(c(-3, 0, 3), c("#2166AC", "white", "#B2182B")),
              top_annotation = col_ha, right_annotation = row_ha,
              cluster_rows = FALSE, cluster_columns = TRUE,
              show_row_names = FALSE, show_column_names = FALSE,
              row_split = factor(ifelse(rownames(Mz) %in% CORE, "core candidates", "background"),
                                 levels = c("core candidates", "background")),
              row_title_gp = gpar(fontsize = 8),
              column_title = sprintf("Knockout impact spectrum: %d knockouts x %d affected genes",
                                     nrow(Mz), ncol(Mz)),
              column_title_gp = gpar(fontsize = 10, fontface = "bold"),
              heatmap_legend_param = list(title_gp = gpar(fontsize = 8),
                                          labels_gp = gpar(fontsize = 7)))
draw(ht, merge_legend = TRUE)
pdf(file.path(fig, "Fig10_knockout_impact_spectrum.pdf"), width = 9.2, height = 6.2)
draw(ht, merge_legend = TRUE); dev.off()
png(file.path(fig, "Fig10_knockout_impact_spectrum.png"), width = 9.2, height = 6.2,
    units = "in", res = 300); draw(ht, merge_legend = TRUE); dev.off()
# 可编辑矢量（ComplexHeatmap 的 svg 走 grid，仍然文本可选可改）
svglite::svglite(file.path(fig, "Fig10_knockout_impact_spectrum.svg"),
                 width = 9.2, height = 6.2)
draw(ht, merge_legend = TRUE); dev.off()

# ---------- Fig11 影响景观（用官方 perturbationMap 的相似度轴） ----------
D <- raw$perturbationDistances
D[is.na(D)] <- 0
S <- raw$perturbationDirections
P <- t(apply(D, 1, rank)) / ncol(D) * as.matrix(S)[, colnames(D), drop = FALSE]
dimnames(P) <- dimnames(D)
selfidx <- match(rownames(P), colnames(P))
P[cbind(which(!is.na(selfidx)), selfidx[!is.na(selfidx)])] <- 0
nP <- sqrt(rowSums(P^2))
P <- P[nP > 0, , drop = FALSE]; nP <- nP[nP > 0]
pca <- prcomp(P, center = TRUE, scale. = FALSE)
xy <- as.data.table(pca$x[, 1:2])
xy[, gene := rownames(P)]
xy[, core := gene %in% CORE]
xy[, epi_effect := imp$mean_abs[match(gene, imp$ko_gene)]]
xy <- merge(xy[, !"epi_effect"], imp[, .(gene = ko_gene, epi_effect = mean_abs)], by = "gene")

p11 <- ggplot(xy, aes(PC1, PC2)) +
  geom_point(data = xy[core == FALSE], colour = "grey75", size = 0.9, alpha = 0.7) +
  geom_point(data = xy[core == TRUE], aes(colour = epi_effect), size = 3.2) +
  ggrepel::geom_text_repel(data = xy[core == TRUE], aes(label = gene, colour = epi_effect),
                           size = 3, fontface = "bold", max.overlaps = 20,
                           box.padding = 0.5, min.segment.length = 0.2) +
  scale_colour_gradientn(colours = c("#2166AC", "#F4A582", "#B2182B"), name = "mean |delta|") +
  labs(title = "Knockout impact landscape",
       subtitle = "each point = one virtual knockout; position from its perturbation profile",
       x = sprintf("PC1 (%.1f%%)", 100 * summary(pca)$importance[2, 1]),
       y = sprintf("PC2 (%.1f%%)", 100 * summary(pca)$importance[2, 2])) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey40"),
        panel.grid.minor = element_blank())
save_all(p11, "Fig11_impact_landscape", 6.4, 5.2)

# ---------- Fig12 上皮 vs 其余 ----------
cmp <- melt(spec[!is.na(abs_distance), .(
  epithel = mean(abs_distance[affected_gene %in% EPI]),
  other   = mean(abs_distance[!affected_gene %in% EPI])
), by = ko_gene], id.vars = "ko_gene", variable.name = "set", value.name = "mean_abs")
cmp[, ko_gene := factor(ko_gene, levels = rev(imp[order(mean_abs)]$ko_gene))]
p12 <- ggplot(cmp, aes(x = ko_gene, y = mean_abs, fill = set)) +
  geom_col(position = position_dodge(0.8), width = 0.8) +
  scale_fill_manual(values = c(epithel = "#C44E52", other = "#B8C4CE"),
                    labels = c(epithel = "epithelial/keratinization genes",
                               other = "all other genes"), name = NULL) +
  labs(x = NULL, y = "mean |delta perturbation distance|",
       title = "Does knocking out any gene hit the epithelial program?",
       subtitle = sprintf("%d virtual knockouts; red bars = effect on the %d-gene epithelial program",
                          nrow(imp), length(EPI))) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
        legend.position = "top", plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 7.5, colour = "grey40"),
        panel.grid.major.x = element_blank())
save_all(p12, "Fig12_epithelial_vs_rest", 9, 4.4)

# ---------- Fig13 第三敲除 ----------
k3 <- ko3[!is_target & !is.na(p.adj)]
mk13 <- function(g) {
  d <- k3[target == g][order(-Z)][1:12]
  d[, gene := factor(gene, levels = rev(gene))]
  d[, sig := p.adj < 0.05]
  d[, on_epi := gene %in% EPI]
  ggplot(d, aes(x = gene, y = Z, fill = on_epi, alpha = sig)) +
    geom_col(width = 0.72) + coord_flip() +
    geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.35) +
    scale_fill_manual(values = c(`TRUE` = "#C44E52", `FALSE` = "#9AA7B4"),
                      labels = c(`TRUE` = "epithelial program", `FALSE` = "other"),
                      name = NULL) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.45),
                       labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "n.s."), name = NULL) +
    labs(x = NULL, y = "differential-regulation Z",
         title = sprintf("%s knockout", g)) +
    theme_bw(base_size = 8.5) +
    theme(plot.title = element_text(face = "bold", size = 9),
          legend.position = "bottom", legend.box = "horizontal",
          legend.text = element_text(size = 6.5),
          legend.key.size = unit(3, "mm"),
          legend.margin = margin(0, 0, 0, 0),
          panel.grid.major.y = element_blank())
}
gs13 <- c("P4hb", "S100a13", "Vim")
gs13 <- intersect(gs13, unique(k3$target))
panels <- lapply(gs13, mk13)
p13 <- wrap_plots(panels, nrow = 1) + plot_annotation(
  title = "Additional knockouts reproduce the same keratinocyte-differentiation top list",
  caption = sprintf(paste0("Pale bars: FDR >= 0.05 - after empirical-null calibration no non-target gene reaches ",
                           "significance, so the ranking is shown instead. Red bars: epithelial/keratinization program. ",
                           "%d-gene network (D8)."), nrow(raw$perturbationDistances)),
  theme = theme(plot.title = element_text(face = "bold", size = 10),
                plot.caption = element_text(size = 6.8, colour = "grey40", hjust = 0)))
save_all(p13, "Fig13_third_knockout", 3.5 * length(panels) + 0.8, 4.2)

cat("figures written:\n"); print(list.files(fig, pattern = "^Fig1[0-3]"))
cat("\nD8-04 done\n")

