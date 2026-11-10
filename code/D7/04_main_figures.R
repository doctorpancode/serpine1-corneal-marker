# ============================================================
# D7-04 主图：9 张论文主图（整合 D2-D7 核心结论）
# 全部从已有 results/*.tsv 拼装，不重跑分析
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})
D7 <- "D:/A生信分析/虚拟敲除/第7天执行结果"
D6 <- "D:/A生信分析/虚拟敲除/第6天执行结果"
D5 <- "D:/A生信分析/虚拟敲除/第5天执行结果"
D4 <- "D:/A生信分析/虚拟敲除/第4天执行结果"
fig <- file.path(D7, "figures", "main")
dir.create(fig, showWarnings = FALSE, recursive = TRUE)
rd <- function(...) fread(file.path(...))

comp <- rd(D6, "results", "D6_celltype_composition_wide.tsv")
au   <- rd(D6, "results", "D6_aucell_summary.tsv")
cc   <- rd(D6, "results", "D6_cellchat_inflam_to_vessel.tsv")
dec  <- rd(D7, "results", "D7_deconvolution_weights.tsv")
kor  <- rd(D7, "results", "D7_knockout_rwr.tsv")
cor2 <- rd(D7, "results", "D7_coregene_vs_celltype_cor.tsv")
enr  <- rd(D7, "results", "D7_enrichment_all.tsv")
gs   <- rd(D5, "results", "D5_group_stats.tsv")
sc   <- rd(D4, "results", "D4_candidate_scorecard.tsv")

# ---------- Fig1: 研究流程（示意，用 ggplot 画 flow）----------
flow <- data.table(
  x = c(1, 2, 3, 4, 5),
  y = c(1, 1, 1, 1, 1),
  step = c("Bulk discovery\n(GSE241715/GSE89537)\n167 candidate genes",
           "WGCNA + stemness\n(D2: 6 modules)",
           "Machine learning\n(D4: LASSO/RF/SVM-RFE)",
           "Single-cell atlas\n(D6: 11 cell types)",
           "Virtual knockout\n(D7: scTenifoldKnk)")
)
p1 <- ggplot(flow, aes(x = x, y = y)) +
  geom_segment(data = flow[1:4], aes(xend = x + 0.72, yend = y),
               arrow = arrow(length = unit(3, "mm")), linewidth = 0.7, colour = "#444444") +
  geom_label(aes(label = step), size = 3, label.padding = unit(0.45, "lines"),
             fill = c("#DEEBF7", "#DEEBF7", "#DEEBF7", "#FEE0D2", "#FEE0D2"),
             label.size = 0.3) +
  coord_cartesian(xlim = c(0.55, 5.55), ylim = c(0.78, 1.22), clip = "off") +
  theme_void(base_size = 11) +
  labs(title = "Study design") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        plot.margin = margin(6, 6, 6, 6))
ggsave(file.path(fig, "Fig1_study_design.png"), p1, width = 14, height = 2.9, dpi = 300)
ggsave(file.path(fig, "Fig1_study_design.pdf"), p1, width = 14, height = 2.9)

# ---------- Fig2: bulk 表达谱（四组小提琴）----------
gs_long <- copy(gs[, .(gene, group, expr = mean)])
gs_long[, gene := factor(gene, levels = gs[group == "Control"][order(-mean)]$gene)]
lvl <- c("Control", "Keratoconus", "Viral", "Bacterial_fungal")
gs_long <- gs_long[group %in% lvl]
p2 <- ggplot(gs_long, aes(x = gene, y = expr, fill = factor(group, levels = lvl))) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  scale_fill_manual(values = c("#4C72B0", "#55A868", "#DD8452", "#C44E52"), name = NULL) +
  coord_flip() +
  labs(x = NULL, y = "Mean logCPM (GSE241715, n = 21)",
       title = "Core genes across corneal disease groups") +
  theme_bw(base_size = 10) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))
ggsave(file.path(fig, "Fig2_bulk_expression.png"), p2, width = 9, height = 6.5, dpi = 300)
ggsave(file.path(fig, "Fig2_bulk_expression.pdf"), p2, width = 9, height = 6.5)

# ---------- Fig3: ML 收敛 ----------
p3 <- ggplot(sc, aes(x = reorder(gene, score), y = score,
                     fill = factor(score))) +
  geom_col(width = 0.72) + coord_flip() +
  scale_fill_brewer(palette = "YlOrRd", guide = "none") +
  geom_text(aes(label = score), hjust = -0.35, size = 3) +
  coord_cartesian(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Evidence score (0-6)", title = "Candidate prioritisation scorecard") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        plot.margin = margin(5, 30, 5, 5))
ggsave(file.path(fig, "Fig3_scorecard.png"), p3, width = 7.5, height = 5, dpi = 300)
ggsave(file.path(fig, "Fig3_scorecard.pdf"), p3, width = 7.5, height = 5)

# ---------- Fig4: 单细胞组成 ----------
setnames(comp, c("celltype", "n_NaOH", "n_uninjury", "prop_NaOH", "prop_uninjury"))
cl <- melt(comp[, .(celltype, NaOH = prop_NaOH, uninjury = prop_uninjury)],
           id.vars = "celltype", variable.name = "sample", value.name = "prop")
cl[, celltype := gsub("_", " ", celltype)]
p4 <- ggplot(cl, aes(x = reorder(celltype, prop), y = prop, fill = sample)) +
  geom_col(position = position_dodge(0.8), width = 0.75) + coord_flip() +
  scale_fill_manual(values = c(NaOH = "#C44E52", uninjury = "#4C72B0"), name = NULL) +
  labs(x = NULL, y = "Fraction of cells", title = "Cell composition shift after alkali burn") +
  theme_bw(base_size = 10) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))
ggsave(file.path(fig, "Fig4_cell_composition.png"), p4, width = 8.5, height = 5.2, dpi = 300)
ggsave(file.path(fig, "Fig4_cell_composition.pdf"), p4, width = 8.5, height = 5.2)

# ---------- Fig5: SERPINE1 在内皮 ----------
a5 <- au[gene == "Serpine1" & n >= 20]
a5[, ct := gsub("_", " ", celltype)]
p5 <- ggplot(a5, aes(x = reorder(ct, mean_auc), y = mean_auc, fill = sample)) +
  geom_col(position = position_dodge(0.8), width = 0.75) + coord_flip() +
  scale_fill_manual(values = c(NaOH = "#C44E52", uninjury = "#4C72B0"), name = NULL) +
  labs(x = NULL, y = "Mean AUCell score", title = "SERPINE1 activity is endothelium-restricted") +
  theme_bw(base_size = 10) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))
ggsave(file.path(fig, "Fig5_serpine1_aucell.png"), p5, width = 8.5, height = 5.2, dpi = 300)
ggsave(file.path(fig, "Fig5_serpine1_aucell.pdf"), p5, width = 8.5, height = 5.2)

# ---------- Fig6: 通讯轴 ----------
cc[, pair := paste0(ligand, "-", receptor)]
t6 <- head(cc[order(-strength)], 12)
t6[, pair := factor(paste0(pair, " (", sub("_", " ", source), "->", sub("_", " ", target), ")"),
                    levels = rev(paste0(pair, " (", sub("_", " ", source), "->", sub("_", " ", target), ")")))]
p6 <- ggplot(t6, aes(x = pair, y = strength, fill = target)) +
  geom_col(width = 0.72) + coord_flip() +
  scale_fill_manual(values = c(Endothelium = "#55A868", Endothelium_Lymphatic = "#DD8452"), name = "target") +
  labs(x = NULL, y = "Communication strength", title = "Inflammation-to-vessel ligand-receptor axes") +
  theme_bw(base_size = 10) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))
ggsave(file.path(fig, "Fig6_cellchat.png"), p6, width = 10, height = 5.4, dpi = 300)
ggsave(file.path(fig, "Fig6_cellchat.pdf"), p6, width = 10, height = 5.4)

# ---------- Fig7: 解卷积 ----------
dec[, celltype := gsub("_", " ", celltype)]
ordw <- dec[, .(m = mean(weight)), by = .(celltype)][order(-m)]$celltype
dec[, celltype := factor(celltype, levels = ordw)]
p7 <- ggplot(dec, aes(x = sample, y = weight, fill = celltype)) +
  geom_col(width = 0.85) +
  labs(x = NULL, y = "Estimated fraction (NNLS)", fill = "cell type",
       title = "Bulk deconvolution with the single-cell signature") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
        legend.position = "right", plot.title = element_text(face = "bold"))
ggsave(file.path(fig, "Fig7_deconvolution.png"), p7, width = 11, height = 5.4, dpi = 300)
ggsave(file.path(fig, "Fig7_deconvolution.pdf"), p7, width = 11, height = 5.4)

# ---------- Fig8: 虚拟敲除 ----------
mk8 <- function(dt, nm, n = 12) {
  d <- dt[target == nm & gene != nm][order(-abs(delta))][1:n]
  d[, gene := factor(gene, levels = rev(gene))]
  ggplot(d, aes(x = gene, y = delta, fill = delta > 0)) +
    geom_col(width = 0.72) + coord_flip() +
    scale_fill_manual(values = c(`TRUE` = "#C44E52", `FALSE` = "#4C72B0"), guide = "none") +
    labs(x = NULL, y = "RWR influence change", title = sprintf("%s knockout", nm)) +
    theme_bw(base_size = 9) +
    theme(plot.title = element_text(face = "bold"))
}
p8 <- mk8(kor, "Serpine1") | mk8(kor, "Tmsb4x")
ggsave(file.path(fig, "Fig8_virtual_knockout.png"), p8, width = 11, height = 4.6, dpi = 300)
ggsave(file.path(fig, "Fig8_virtual_knockout.pdf"), p8, width = 11, height = 4.6)

# ---------- Fig9: 机制模型 ----------
nodes <- data.table(
  x = c(1, 2.2, 3.4, 4.6, 5.8),
  y = 3,
  lab = c("Alkali burn /\nsuture injury",
          "Neutrophil &\nmacrophage influx",
          "SERPINE1\n(endothelium)",
          "TMSB4X\n(epithelium)",
          "Neovascularization\n& lymphangiogenesis")
)
notes <- data.table(
  x = c(3.4, 4.6),
  y = c(2.25, 1.85),
  lab = c("Lgals3-Cd44 / Spp1-Cd44 / Il1b-Il1r1",
          "ECM organisation / PI3K-Akt"),
  col = c("#C44E52", "#4C72B0")
)
p9 <- ggplot() +
  geom_segment(data = nodes[1:4], aes(x = x, xend = x + 0.82, y = y, yend = y),
               arrow = arrow(length = unit(3, "mm")), linewidth = 0.8, colour = "#555555") +
  geom_label(data = nodes, aes(x = x, y = y, label = lab), size = 3.1,
             label.padding = unit(0.45, "lines"), label.size = 0.3,
             fill = c("#FEE0D2", "#FFF7BC", "#DEEBF7", "#DEEBF7", "#E5F5E0")) +
  geom_text(data = notes, aes(x = x, y = y, label = lab, colour = col), size = 3.6) +
  scale_colour_identity() +
  coord_cartesian(xlim = c(0.45, 6.55), ylim = c(1.55, 3.45), clip = "off") +
  theme_void(base_size = 11) +
  labs(title = "Proposed SERPINE1-TMSB4X axis in corneal neovascularization") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
        plot.margin = margin(10, 10, 10, 10))
ggsave(file.path(fig, "Fig9_mechanism_model.png"), p9, width = 14, height = 4.0, dpi = 300)
ggsave(file.path(fig, "Fig9_mechanism_model.pdf"), p9, width = 14, height = 4.0)

cat("main figures written to", fig, "\n")
print(list.files(fig))
cat("\nD7-04 done\n")