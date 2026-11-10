suppressPackageStartupMessages({library(Seurat)})
obj <- readRDS("第6天执行结果/data/D6_seurat_GSE275678_final.rds")
cat("dim:", dim(obj), "\n")
m <- obj@meta.data
cat("meta cols:", paste(colnames(m), collapse=" | "), "\n")
um <- Embeddings(obj, "umap")
cat("umap dim:", dim(um), "\n")
out <- data.frame(cell = rownames(m), UMAP_1 = um[,1], UMAP_2 = um[,2],
                  m, check.names = FALSE)
write.table(out, "文章/figs_v2/out/D6_umap_meta.tsv", sep="\t",
            quote = FALSE, row.names = FALSE)
cat("written\n")
print(table(m[["celltype"]]))