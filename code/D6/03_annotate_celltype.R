# ============================================================
# D6-03 细胞类型注释（规则式：marker 平均表达 z-score 最大者）
# 并做 cluster 的 top marker 佐证
# ============================================================
suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(ggplot2); library(patchwork)
})
set.seed(20261109)
out <- "D:/A生信分析/虚拟敲除/第6天执行结果"
fig <- file.path(out, "figures")

o <- readRDS(file.path(out, "data", "D6_seurat_GSE275678_annot_input.rds"))
sc <- fread(file.path(out, "results", "D6_cluster_marker_scores.tsv"))
score_cols <- grep("1$", names(sc), value = TRUE)
score_cols <- setdiff(score_cols, c("n", "cluster"))

# z-score 跨 cluster 标准化后取最大
z <- as.matrix(sc[, ..score_cols])
z <- scale(z)
colnames(z) <- sub("1$", "", colnames(z))
lab <- colnames(z)[apply(z, 1, which.max)]
sc[, auto_label := lab]
sc[, margin := apply(z, 1, function(v) sort(v, decreasing = TRUE)[1] - sort(v, decreasing = TRUE)[2])]

# ---- 每个 cluster 的 top marker（用 Seurat 内置检验）----
o <- JoinLayers(o)   # Seurat v5 多 layer 必须先合并，否则 FindAllMarkers 不跑检验
mk <- as.data.table(FindAllMarkers(o, only.pos = TRUE, min.pct = 0.25,
                                   logfc.threshold = 0.5, verbose = FALSE))
fwrite(mk, file.path(out, "results", "D6_cluster_top_markers.tsv"), sep = "\t")
top5 <- mk[order(cluster, -avg_log2FC)][, head(.SD, 5), by = cluster]
fwrite(top5, file.path(out, "results", "D6_cluster_top5_markers.tsv"), sep = "\t")
print(top5[, .(cluster, gene, avg_log2FC = round(avg_log2FC, 2), pct.1)])

fwrite(sc, file.path(out, "results", "D6_cluster_annotation_auto.tsv"), sep = "\t")

# auto_label 不可靠（z-score 跨 cluster 标准化后被单一 marker 主导），此处只保存对象供 03b 手工注释
saveRDS(o, file.path(out, "data", "D6_seurat_GSE275678_annotated.rds"))
cat("\nD6-03 done (marker only)\n")
