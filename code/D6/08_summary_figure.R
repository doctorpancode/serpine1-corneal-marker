# ============================================================
# D6-08 汇总图：单细胞五件套的核心结论（论文 Fig5 候选）
# 输入全部来自本目录 results/*.tsv，不重跑 Seurat
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
out <- "D:/A生信分析/虚拟敲除/第6天执行结果"
R <- function(f) fread(file.path(out, "results", f))
fig <- file.path(out, "figures")

# ---- A: 细胞组成变化（NaOH / uninjury 比例 + log2FC）----
comp <- R("D6_celltype_composition_wide.tsv")
setnames(comp, c("celltype", "n_NaOH", "n_uninjury", "prop_NaOH", "prop_uninjury"))
comp[, `:=`(label = gsub("_", " ", celltype),
            ratio = (prop_NaOH + 2e-4) / (prop_uninjury + 2e-4))]
comp[, log2ratio := log2(ratio)]
comp <- comp[order(-prop_NaOH)]
comp[, label := factor(label, levels = rev(label))]

pA <- ggplot(comp, aes(x = label, y = log2ratio, fill = log2ratio > 0)) +
  geom_col(width = 0.72) + coord_flip() + geom_hline(yintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = c(`TRUE` = "#C44E52", `FALSE` = "#4C72B0"), guide = "none") +
  labs(x = NULL, y = "log2 (lesion / normal) of cell fraction",
       title = "A  Cell composition shift (NaOH vs uninjury)") +
  theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold", hjust = 0))

# ---- B: Serpine1 AUCell 活性（内皮最强）----
au <- R("D6_aucell_summary.tsv")
au1 <- au[gene == "Serpine1" & n >= 20]
au1[, label := gsub("_", " ", celltype)]
au1[, label := factor(label, levels = rev(unique(label[order(-mean_auc)])))]
pB <- ggplot(au1, aes(x = label, y = mean_auc, fill = sample)) +
  geom_col(position = position_dodge(0.78), width = 0.72) +
  scale_fill_manual(values = c(NaOH = "#C44E52", uninjury = "#4C72B0"), name = NULL) +
  coord_flip() +
  labs(x = NULL, y = "Mean AUCell score (Serpine1 regulon)",
       title = "B  Serpine1 activity by cell type") +
  theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold", hjust = 0),
                                   legend.position = "top")

# ---- C: 炎症→血管通讯（top 12）----
cc <- R("D6_cellchat_inflam_to_vessel.tsv")
cc[, pair := paste0(ligand, "-", receptor, "  ", gsub("_", " ", source), " -> ", gsub("_", " ", target))]
top <- head(cc[order(-strength)], 12)
top[, pair := factor(pair, levels = rev(pair))]
pC <- ggplot(top, aes(x = pair, y = strength, fill = target)) +
  geom_col(width = 0.72) + coord_flip() +
  scale_fill_manual(values = c(Endothelium_Lymphatic = "#DD8452", Endothelium = "#55A868"), name = "target") +
  labs(x = NULL, y = "Communication strength", title = "C  Inflammation-to-vessel axis") +
  theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold", hjust = 0),
                                   legend.position = "top")

# ---- D: 干性丢失（NaOH vs uninjury）----
st <- R("D6_stemness_by_celltype.tsv")
st[, label := gsub("_", " ", celltype)]
pD <- ggplot(st[!is.na(mean_stem)], aes(x = label, y = mean_stem, fill = sample)) +
  geom_col(position = position_dodge(0.78), width = 0.72) +
  scale_fill_manual(values = c(NaOH = "#C44E52", uninjury = "#4C72B0"), name = NULL) +
  coord_flip() +
  labs(x = NULL, y = "Mean stemness module score",
       title = "D  Stemness program by cell type") +
  theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold", hjust = 0),
                                   legend.position = "top")

p <- (pA | pB) / (pC | pD) +
  plot_annotation(title = "Single-cell atlas of corneal neovascularization (GSE275678)",
                  subtitle = "19,493 genes x 15,873 cells; 11 cell types; NaOH burn vs uninjured cornea",
                  theme = theme(plot.title = element_text(face = "bold", size = 14)))
ggsave(file.path(fig, "D6_Fig12_summary_panel.png"), p, width = 15, height = 11, dpi = 300)
ggsave(file.path(fig, "D6_Fig12_summary_panel.pdf"), p, width = 15, height = 11)

fwrite(comp, file.path(out, "results", "D6_summary_composition.tsv"), sep = "\t")
cat("\nD6-08 done\n")