# ============================================================
# D6-03b 用 top marker 做人工/半自动注释（替代 auto_label，后者不可靠）
# ============================================================
suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(ggplot2)
})
out <- "D:/A生信分析/虚拟敲除/第6天执行结果"
fig <- file.path(out, "figures")

o <- readRDS(file.path(out, "data", "D6_seurat_GSE275678_annotated.rds"))
top5 <- fread(file.path(out, "results", "D6_cluster_top5_markers.tsv"))
sc   <- fread(file.path(out, "results", "D6_cluster_annotation_auto.tsv"))

cat("\n===== 每个 cluster 的 top5 marker =====\n")
for (k in sort(unique(top5$cluster))) {
  cat(sprintf("C%-2s n=%-5s | %s\n", k, sc[cluster == k]$n,
              paste(top5[cluster == k]$gene, collapse = ", ")))
}

# ---- 人工注释表（依据 top5 marker + 经典角膜细胞类型）----
ann <- c(
  "0"  = "Epithelium_Basal",     # Aldh3a1/Nppb 干/基底
  "1"  = "Epithelium_Wing",
  "2"  = "Epithelium_Basal",
  "3"  = "Fibroblast",
  "4"  = "Neutrophil",
  "5"  = "Epithelium_Superficial",
  "6"  = "Macrophage",
  "7"  = "Epithelium_Basal",
  "8"  = "Fibroblast",
  "9"  = "Epithelium_Wing",
  "10" = "Tcell_NK",
  "11" = "Epithelium_Superficial",
  "12" = "Epithelium_Basal",
  "13" = "Fibroblast",
  "14" = "Epithelium_Wing",
  "15" = "Macrophage",
  "16" = "Endothelium",
  "17" = "Pericyte_SMC",
  "18" = "Epithelium_Basal",
  "19" = "Epithelium_Wing",
  "20" = "Fibroblast"
)
o$celltype <- ifelse(as.character(o$seurat_clusters) %in% names(ann),
                     ann[as.character(o$seurat_clusters)], "Unknown")
cat("\n===== celltype counts =====\n"); print(table(o$celltype))

md <- as.data.table(o@meta.data)
prop <- md[, .N, by = .(sample, celltype)][, frac := N / sum(N), by = sample]
fwrite(prop, file.path(out, "results", "D6_celltype_proportion_by_sample.tsv"), sep = "\t")
print(dcast(prop, celltype ~ sample, value.var = "frac"))

saveRDS(o, file.path(out, "data", "D6_seurat_GSE275678_final.rds"))

p1 <- DimPlot(o, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) +
  ggtitle("Cell types (GSE275678)") + theme(legend.position = "right")
ggsave(file.path(fig, "D6_Fig3_UMAP_celltype.png"), p1, width = 10.5, height = 6.5, dpi = 300)
ggsave(file.path(fig, "D6_Fig3_UMAP_celltype.pdf"), p1, width = 10.5, height = 6.5)

p2 <- ggplot(prop, aes(x = reorder(celltype, -frac), y = frac, fill = sample)) +
  geom_col(position = "dodge") + coord_flip() +
  labs(x = "", y = "Fraction of cells",
       title = "Cell composition: uninjured vs alkali-burn cornea") +
  scale_fill_manual(values = c(uninjury = "#4C72B0", NaOH = "#C44E52")) +
  theme_bw(base_size = 12)
ggsave(file.path(fig, "D6_Fig4_celltype_composition.png"), p2, width = 9.5, height = 5.5, dpi = 300)
ggsave(file.path(fig, "D6_Fig4_celltype_composition.pdf"), p2, width = 9.5, height = 5.5)
cat("\nD6-03b done\n")