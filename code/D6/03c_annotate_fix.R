# ============================================================
# D6-03c 修正注释：基于真实 top5 marker 重新判定
# 关键发现：C3=中性粒细胞(S100a8/Cxcr2/Stfa2l1)、C14=单核/巨噬(Msr1/Ccr2/Ms4a6)、
#          C6/C8/C12=增殖(SMC2/Cenpf/Nusap1/Mcm)、C18=周细胞(Rgs5/Rgs4/Vtn)、
#          C16=血管内皮(Sele/Ptprb/Sox18/Myct1)、C17=施万(Plp1/Mpz/Gfra3)、
#          C2/C11=上皮(Nppb/Aldh3a1/Krt19/Ereg)、C13=基质(Kera/Angptl7/Pcsk2)
# 保留：C0/C5/C7/C9/C19 = 上皮亚群
# ============================================================
suppressPackageStartupMessages({ library(Seurat); library(data.table); library(ggplot2) })
out <- "D:/A生信分析/虚拟敲除/第6天执行结果"; fig <- file.path(out, "figures")

o <- readRDS(file.path(out, "data", "D6_seurat_GSE275678_annotated.rds"))
top5 <- fread(file.path(out, "results", "D6_cluster_top5_markers.tsv"))
sc <- fread(file.path(out, "results", "D6_cluster_annotation_auto.tsv"))

ann <- c(
  "0"  = "Epithelium_Basal",      # Nppb/Necab1/Aldh3a1
  "1"  = "Epithelium_Wing",       # Hpgd/Gpha2/Car12/Ereg
  "2"  = "Epithelium_Basal",      # Cyp2a5/Sprr2f/Krt19
  "3"  = "Neutrophil",            # S100a8/Cxcr2/Stfa2l1/Acod1
  "4"  = "Epithelium_Proliferating", # Hist1 组蛋白簇
  "5"  = "Epithelium_Wing",       # Ly6g6c/Cldn17/Muc15
  "6"  = "Epithelium_Proliferating", # Mns1/Cenpf/Nusap1/Kif2c/Cdkn3
  "7"  = "Epithelium_Basal",      # Prg4/Rspo1/Dio2/Ecrg4
  "8"  = "Epithelium_Proliferating", # Ung/Hells/Chek1/Mcm3/Chaf1b
  "9"  = "Epithelium_Superficial",# Dsg1a/Dsg3/Mal
  "10" = "Endothelium_Lymphatic", # Col4a4/Prss12/Wfdc18
  "11" = "Epithelium_Wing",       # Serpinb3a/Esp6/Car12
  "12" = "Epithelium_Proliferating", # Dscc1/Mcm2/Mcm5/Lig1
  "13" = "Fibroblast",            # Kera/Angptl7/Pcsk2/Tacr3
  "14" = "Macrophage",            # Msr1/Ccr2/Ms4a6/Ms4a6c
  "15" = "Epithelium_Basal",      # Atp1b1/Wnt6/Id1/Tppp3
  "16" = "Endothelium",           # Sele/Ptprb/Sox18/Myct1
  "17" = "Schwann",               # Cdh19/Mpz/Gfra3/Gpr37l1/Kcna1
  "18" = "Pericyte_SMC",          # Rgs5/Rgs4/Vtn/Higd1b
  "19" = "Epithelium_Superficial",# Lce3a/Clca4b
  "20" = "Epithelium_Wing"        # Thsd4/Hpgd/Ccdc3 (C1 的小卫星)
)
o$celltype <- ifelse(as.character(o$seurat_clusters) %in% names(ann),
                     ann[as.character(o$seurat_clusters)], "Unknown")
cat("\n===== celltype counts =====\n"); print(table(o$celltype))

md <- as.data.table(o@meta.data)
prop <- md[, .N, by = .(sample, celltype)][, frac := N / sum(N), by = sample]
fwrite(prop, file.path(out, "results", "D6_celltype_proportion_by_sample.tsv"), sep = "\t")
w <- dcast(prop, celltype ~ sample, value.var = c("N", "frac"), fill = 0)
fwrite(w, file.path(out, "results", "D6_celltype_composition_wide.tsv"), sep = "\t")
print(dcast(prop, celltype ~ sample, value.var = "frac", fill = 0))

saveRDS(o, file.path(out, "data", "D6_seurat_GSE275678_final.rds"))

p1 <- DimPlot(o, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) +
  ggtitle("Cell types (GSE275678, alkali burn vs uninjured)") +
  theme(legend.text = element_text(size = 9))
ggsave(file.path(fig, "D6_Fig3_UMAP_celltype.png"), p1, width = 11.5, height = 6.5, dpi = 300)
ggsave(file.path(fig, "D6_Fig3_UMAP_celltype.pdf"), p1, width = 11.5, height = 6.5)

prop[, celltype := factor(celltype, levels = prop[sample == "NaOH"][order(-frac)]$celltype)]
p2 <- ggplot(prop, aes(x = celltype, y = frac, fill = sample)) +
  geom_col(position = "dodge") + coord_flip() +
  labs(x = "", y = "Fraction of cells",
       title = "Cell composition shift after alkali burn") +
  scale_fill_manual(values = c(uninjury = "#4C72B0", NaOH = "#C44E52")) +
  theme_bw(base_size = 12)
ggsave(file.path(fig, "D6_Fig4_celltype_composition.png"), p2, width = 9.5, height = 5.8, dpi = 300)
ggsave(file.path(fig, "D6_Fig4_celltype_composition.pdf"), p2, width = 9.5, height = 5.8)
cat("\nD6-03c done\n")