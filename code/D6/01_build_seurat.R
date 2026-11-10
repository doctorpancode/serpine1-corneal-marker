# ============================================================
# D6-01 单细胞图谱构建（Seurat）
# 数据：GSE275678（鼠 碱烧伤 vs 未损伤，2 样本）为主
#        GSE294213（鼠 缝线 CNV，3 样本）为验证
# 输出：D:/A生信分析/虚拟敲除/第6天执行结果/
# ============================================================
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table); library(ggplot2)
})
set.seed(20261109)
out <- "D:/A生信分析/虚拟敲除/第6天执行结果"
raw <- "D:/A生信分析/虚拟敲除/D1执行结果/D1_CNV_角膜新生血管_数据分析/data/raw"
dir.create(file.path(out, "data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "results"), recursive = TRUE, showWarnings = FALSE)

read10x_one <- function(gse, sample_txt) {
  d <- file.path(raw, gse, "extracted"); fs <- list.files(d, pattern = sample_txt, full.names = TRUE)
  m <- ReadMtx(
    mtx   = fs[grepl("matrix", fs)][1],
    cells = fs[grepl("barcodes", fs)][1],
    features = fs[grepl("features", fs)][1],
    feature.column = 2
  )
  m
}

build_one <- function(gse, samples, prefix) {
  cat("\n==========", gse, "==========\n")
  objs <- list()
  for (s in samples) {
    m <- read10x_one(gse, s)
    cat(s, dim(m), "\n")
    o <- CreateSeuratObject(m, project = paste0(prefix, "_", s), min.cells = 3, min.features = 200)
    o$sample <- s; o$gse <- gse
    objs[[s]] <- o
  }
  o <- merge(objs[[1]], objs[-1], add.cell.ids = samples)
  o[["percent.mt"]] <- PercentageFeatureSet(o, pattern = "^mt-")
  cat("merged:", dim(o), "\n")
  o
}

# ---- GSE275678：uninjury vs NaOH（2 样本，主分析）----
s275678 <- build_one("GSE275678", c("uninjury", "NaOH"), "GSE275678")
s275678 <- subset(s275678, subset = nFeature_RNA > 300 & nFeature_RNA < 8000 & percent.mt < 15)
cat("GSE275678 after QC:", dim(s275678), "\n")

# ---- GSE294213：nosut / sut / dultreatment（3 样本，验证）----
s294213 <- build_one("GSE294213", c("nosut", "sut", "dultreatment"), "GSE294213")
s294213 <- subset(s294213, subset = nFeature_RNA > 300 & nFeature_RNA < 8000 & percent.mt < 15)
cat("GSE294213 after QC:", dim(s294213), "\n")

saveRDS(s275678, file.path(out, "data", "D6_seurat_GSE275678_qc.rds"))
saveRDS(s294213, file.path(out, "data", "D6_seurat_GSE294213_qc.rds"))

qc_fun <- function(o, g) {
  md <- as.data.table(o@meta.data)
  md[, .(gse = g, n_cells = .N,
         med_features = as.numeric(median(nFeature_RNA)),
         med_counts   = as.numeric(median(nCount_RNA)),
         med_mt       = as.numeric(median(percent.mt))), by = sample]
}
qc <- rbind(qc_fun(s275678, "GSE275678"), qc_fun(s294213, "GSE294213"))
fwrite(qc, file.path(out, "results", "D6_qc_after_filter.tsv"), sep = "\t")
print(qc)
cat("\nD6-01 done\n")
