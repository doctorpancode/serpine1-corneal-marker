# ============================================================
# D6-05 干性评分（复用 D2 的 61 基因集）+ 拟时序（自实现）
# ============================================================
suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(ggplot2); library(patchwork)
  library(igraph)
})
set.seed(20261109)
out <- "D:/A生信分析/虚拟敲除/第6天执行结果"; fig <- file.path(out, "figures")

o <- readRDS(file.path(out, "data", "D6_seurat_GSE275678_final.rds"))

# ---- 1) 干性评分：D2 的 61 基因集，鼠符号映射（首字母大写）----
gs <- fread("D:/A生信分析/虚拟敲除/第2天执行结果/data/stemness_genesets.tsv")
cat("D2 gene sets:", nrow(gs), "unique:", length(unique(gs$gene)), "\n")
mouse <- function(x) paste0(toupper(substr(x,1,1)), tolower(substring(x,2)))
gs[, mouse_gene := mouse(gene)]
gsm <- unique(gs$mouse_gene)
present <- intersect(gsm, rownames(o))
cat("mapped to mouse and present:", length(present), "/", length(gsm), "\n")
miss <- setdiff(gsm, rownames(o))
cat("missing:", paste(head(miss, 20), collapse=", "), "\n")

o <- AddModuleScore(o, features = list(present), name = "Stemness", seed = 20261109)
md <- as.data.table(o@meta.data)
md[, celltype := o$celltype]
st <- md[, .(mean_stem = mean(Stemness1), med_stem = median(Stemness1), n = .N), by = .(celltype, sample)]
fwrite(st, file.path(out, "results", "D6_stemness_by_celltype.tsv"), sep = "\t")
print(dcast(st, celltype ~ sample, value.var = "mean_stem"))

p_s <- ggplot(st, aes(x = reorder(celltype, mean_stem), y = mean_stem, fill = sample)) +
  geom_col(position = "dodge") + coord_flip() +
  scale_fill_manual(values = c(uninjury = "#4C72B0", NaOH = "#C44E52")) +
  labs(x = "", y = "Stemness score (D2 61-gene set)", title = "Stemness across cell types") +
  theme_bw(base_size = 12)
ggsave(file.path(fig, "D6_Fig8_stemness_celltype.png"), p_s, width = 9.5, height = 5.8, dpi = 300)
ggsave(file.path(fig, "D6_Fig8_stemness_celltype.pdf"), p_s, width = 9.5, height = 5.8)

p_v <- VlnPlot(o, features = "Stemness1", group.by = "celltype", pt.size = 0) +
  labs(x = "", y = "Stemness score", title = "Stemness distribution by cell type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig, "D6_Fig8B_stemness_violin.png"), p_v, width = 10, height = 5.5, dpi = 300)
ggsave(file.path(fig, "D6_Fig8B_stemness_violin.pdf"), p_v, width = 10, height = 5.5)

saveRDS(o, file.path(out, "data", "D6_seurat_GSE275678_stem.rds"))

# ---- 2) 拟时序（自实现：以 basal 上皮为根的扩散伪时间）----
# 只用上皮亚群，避免跨谱系误导
epi <- subset(o, subset = celltype %in% c("Epithelium_Basal","Epithelium_Wing",
                                          "Epithelium_Superficial","Epithelium_Proliferating"))
cat("\nepithelial cells:", ncol(epi), "\n")
epi <- NormalizeData(epi, verbose = FALSE)
epi <- FindVariableFeatures(epi, nfeatures = 1500, verbose = FALSE)
epi <- ScaleData(epi, verbose = FALSE)
epi <- RunPCA(epi, npcs = 20, verbose = FALSE)
epi <- FindNeighbors(epi, dims = 1:20, verbose = FALSE)
epi <- RunUMAP(epi, dims = 1:20, verbose = FALSE)

# 根细胞：basal 中干性最高的 5%
emb <- Embeddings(epi, "pca")[, 1:20]
root_cells <- as.data.table(epi@meta.data)[, .(cell = rownames(.SD), Stemness1, celltype)]
root_cells <- root_cells[celltype == "Epithelium_Basal"][order(-Stemness1)]
roots <- head(root_cells$cell, max(50, round(0.05 * ncol(epi))))

# ---- 伪时间：UMAP 2D 上的 kNN 图 Dijkstra 扩散 ----
# 注意几个坑：
#  1) Seurat FindNeighbors 返回 dgCMatrix（值 1 = 是否近邻），graph_from_adj_list 会退化 → 必须显式算 kNN
#  2) FNN::get.knn 对 k > n-1 会报错，且 nn.index 里可能出现自环
#  3) distances() 的起点用顶点名（cell barcode），不是索引
#  4) 若图不连通，distances 会返回 Inf → apply(...,min) 得 Inf，rescale 后全 NaN
um <- Embeddings(epi, "umap")[, 1:2]
rownames(um) <- colnames(epi)
n <- nrow(um); kN <- min(15, n - 1)
dd_all <- as.matrix(dist(um))
el <- rbindlist(lapply(seq_len(n), function(i) {
  o <- order(dd_all[i, ])[2:(kN + 1)]
  data.table(from = rownames(um)[i], to = rownames(um)[o], w = dd_all[i, o])
}))
g <- graph_from_data_frame(el, directed = FALSE, vertices = data.frame(name = rownames(um)))
g <- simplify(g, edge.attr.comb = "min")

# 保证连通：取最大连通分量，其余细胞丢弃（并在日志里报告数量）
comps <- components(g)
keep_v <- names(comps$membership)[comps$membership == which.max(comps$csize)]
cat("graph vertices:", vcount(g), " components:", comps$no,
    " largest:", length(keep_v), " empty:", sum(comps$csize == 0), "\n")
g <- induced_subgraph(g, keep_v)

roots_use <- intersect(roots, keep_v)
if (!length(roots_use)) roots_use <- head(keep_v, 50)
cat("root cells used:", length(roots_use), "\n")
d <- distances(g, v = roots_use, mode = "all")
pt <- apply(d, 2, min)
cat("finite pt:", sum(is.finite(pt)), "/", length(pt), "\n")
pt[!is.finite(pt)] <- max(pt[is.finite(pt)])
epi$pseudotime <- NA_real_
epi$pseudotime[colnames(epi) %in% names(pt)] <- scales::rescale(pt, to = c(0, 1))
cat("pseudotime range:", range(epi$pseudotime), "\n")

# 与干性、细胞类型的相关性
md2 <- as.data.table(epi@meta.data)[!is.na(pseudotime)]
st_cor <- md2[, .(n = .N, median_pt = median(pseudotime), mean_pt = mean(pseudotime)), by = celltype]
fwrite(st_cor, file.path(out, "results", "D6_pseudotime_by_celltype.tsv"), sep = "\t")
print(st_cor)
ct_cor <- cor.test(md2$pseudotime, md2$Stemness1, method = "spearman")
cat("pseudotime vs stemness: rho =", round(ct_cor$estimate, 3), "P =", signif(ct_cor$p.value, 3), "\n")
fwrite(data.table(comparison = "pseudotime_vs_stemness", rho = ct_cor$estimate,
                  P = ct_cor$p.value, n = nrow(md2)),
       file.path(out, "results", "D6_pseudotime_stemness_cor.tsv"), sep = "\t")

# 主干基因随伪时间
expr <- GetAssayData(epi, layer = "data")
ptdt <- rbindlist(lapply(c("Serpine1","Tmsb4x","Krt3","Krt12","Krt5","Krt14","Krt15"), function(g) {
  if (!g %in% rownames(expr)) return(NULL)
  data.table(gene = g, pseudotime = epi$pseudotime, expr = as.numeric(expr[g, ]))
}))
fwrite(ptdt, file.path(out, "results", "D6_pseudotime_gene_expr.tsv"), sep = "\t")

p_pt <- ggplot(md2, aes(x = pseudotime, fill = celltype)) +
  geom_density(alpha = 0.55) +
  labs(x = "Diffusion pseudotime", y = "Density", title = "Epithelial differentiation trajectory") +
  theme_bw(base_size = 12)
p_um <- FeaturePlot(epi, features = "pseudotime", pt.size = 0.4) +
  scale_colour_gradientn(colours = c("#4575B4","#FFFFBF","#D73027")) + ggtitle("Pseudotime")
ggsave(file.path(fig, "D6_Fig9_pseudotime.png"), wrap_plots(p_um, p_pt, ncol = 2, widths = c(1.1, 1)),
       width = 13, height = 5.2, dpi = 300)
ggsave(file.path(fig, "D6_Fig9_pseudotime.pdf"), wrap_plots(p_um, p_pt, ncol = 2, widths = c(1.1, 1)),
       width = 13, height = 5.2)

saveRDS(epi, file.path(out, "data", "D6_seurat_epi_pseudotime.rds"))
cat("\nD6-05 done\n")