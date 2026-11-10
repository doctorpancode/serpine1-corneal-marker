# ============================================================
# D7-01 反投影解卷积：用 D6 的 11 个细胞类型签名解 GSE241715 bulk
# 方法：NNLS（非负最小二乘），签名 = 各细胞类型平均 log-normalized 表达
# 输出：每个 bulk 样本的 11 个细胞类型比例 + 双主干基因与细胞比例的关联
# ============================================================
suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(ggplot2)
})
set.seed(20261109)
out   <- "D:/A生信分析/虚拟敲除/第7天执行结果"
d6    <- "D:/A生信分析/虚拟敲除/第6天执行结果"
fig   <- file.path(out, "figures")
R7 <- function(f) fread(file.path(out, "results", f))
W1 <- "D:/A生信分析/虚拟敲除/D1执行结果/D1_CNV_角膜新生血管_数据分析/data/processed"

sc <- readRDS(file.path(d6, "data", "D6_seurat_GSE275678_final.rds"))
sc <- JoinLayers(sc)
cat("scRNA:", ncol(sc), "cells,", nrow(sc), "genes\n")
print(table(sc$celltype))

# ---- 1. 小鼠 -> 人 符号映射（首字母大写；人类全大写匹配）----
mg <- rownames(sc)
hg <- toupper(mg)
# Seurat 里小鼠符号已经是首字母大写，toupper 后即人类风格；去掉明显的小鼠特有后缀
hg <- sub("(RIK|Rik)$", "", hg, ignore.case = FALSE)
ok <- !duplicated(hg) & nchar(hg) > 1
cat("genes mapped:", sum(ok), "/", length(mg), "\n")

# ---- 2. 构建签名矩阵（细胞类型 × 基因）----
expr <- GetAssayData(sc, layer = "data")
cts  <- sort(unique(sc$celltype))
sig_list <- lapply(cts, function(ct) {
  cells <- colnames(sc)[sc$celltype == ct]
  if (length(cells) > 600) cells <- sample(cells, 600)
  Matrix::rowMeans(expr[, cells, drop = FALSE])
})
sig <- do.call(cbind, sig_list)
colnames(sig) <- cts
rownames(sig) <- hg
sig <- sig[ok, , drop = FALSE]
sig <- sig[!duplicated(rownames(sig)), , drop = FALSE]
cat("signature:", nrow(sig), "genes x", ncol(sig), "cell types\n")
saveRDS(sig, file.path(out, "data", "D7_signature_matrix.rds"))

# ---- 3. 读 bulk，取交集基因 ----
suppressPackageStartupMessages(library(Biobase))
es  <- readRDS(file.path(W1, "GSE241715_expression.rds"))
bulk <- exprs(es); pd <- as.data.table(pData(es))
cat("bulk:", nrow(bulk), "genes x", ncol(bulk), "samples; range", range(bulk), "\n")
common <- intersect(rownames(sig), rownames(bulk))
cat("common genes:", length(common), "\n")

S <- sig[common, , drop = FALSE]
B <- bulk[common, , drop = FALSE]
# bulk 是 log2CPM，scRNA 是 log-normalized，量纲不同 -> 各自按基因做 z-score 后再比
Sz <- t(scale(t(S))); Sz[is.na(Sz)] <- 0
Bz <- t(scale(t(B))); Bz[is.na(Bz)] <- 0

# ---- 4. NNLS 解卷积 ----
res <- vector("list", ncol(Bz))
for (i in seq_len(ncol(Bz))) {
  y <- as.numeric(Bz[, i])
  A <- Sz                          # genes x celltypes
  D <- crossprod(A); dd <- as.numeric(crossprod(A, y))
  D <- D + diag(1e-8, ncol(D))     # 数值稳定
  qp <- quadprog::solve.QP(D, as.numeric(dd), diag(ncol(D)), rep(0, ncol(D)))
  w <- as.numeric(qp$solution); w[w < 0] <- 0
  if (sum(w) > 0) w <- w / sum(w)
  res[[i]] <- data.table(sample = colnames(Bz)[i], celltype = colnames(Sz), weight = w)
}
dec <- rbindlist(res)
dec <- merge(dec, pd[, .(sample, group)], by = "sample")
fwrite(dec, file.path(out, "results", "D7_deconvolution_weights.tsv"), sep = "\t")
cat("\nNNLS done. weights per sample (first 6):\n")
print(dcast(dec[order(-weight)][, head(.SD, 4), by = sample][order(sample)],
            sample ~ celltype, value.var = "weight"))

# ---- 5. 各组细胞类型比例 ----
grp <- dec[, .(mean_w = mean(weight), sd_w = sd(weight), n = .N), by = .(group, celltype)]
grp <- grp[order(celltype, -mean_w)]
fwrite(grp, file.path(out, "results", "D7_deconvolution_by_group.tsv"), sep = "\t")
print(dcast(grp, celltype ~ group, value.var = "mean_w")[order(-`Bacterial_fungal`)])

# ---- 6. 双主干基因与细胞比例的关联（解卷积回答"表达由谁驱动"）----
core <- c("SERPINE1", "TMSB4X", "P4HB", "S100A13", "SOD1", "CAPG", "VDAC1", "PRDX2", "GJA1", "XIAP")
core <- intersect(core, rownames(B))
wide <- dcast(dec, sample ~ celltype, value.var = "weight")
expr_core <- as.data.table(t(B[core, , drop = FALSE])); setnames(expr_core, core)
expr_core[, sample := colnames(B)]
m <- merge(wide, expr_core, by = "sample")
ct_cols <- setdiff(names(wide), "sample")
cors <- rbindlist(lapply(core, function(g) {
  rbindlist(lapply(ct_cols, function(ct) {
    if (sd(m[[ct]]) == 0 || sd(m[[g]]) == 0) return(NULL)
    ct_ <- suppressWarnings(cor.test(m[[ct]], m[[g]], method = "spearman"))
    data.table(gene = g, celltype = ct, rho = ct_$estimate, P = ct_$p.value, n = nrow(m))
  }))
}))
cors[, FDR := p.adjust(P, "BH")]
cors[, arho := abs(rho)]
setorder(cors, gene, -arho)
fwrite(cors, file.path(out, "results", "D7_coregene_vs_celltype_cor.tsv"), sep = "\t")
cat("\nTop core-gene / cell-type associations (|rho|>0.5, FDR<0.05):\n")
print(cors[abs(rho) > 0.5 & FDR < 0.05][order(gene, -abs(rho))][, .(gene, celltype, rho = round(rho, 3), FDR = signif(FDR, 3))])

# ---- 7. 图 ----
dec[, celltype := factor(celltype)]
p1 <- ggplot(dec, aes(x = sample, y = weight, fill = celltype)) +
  geom_col(width = 0.82) +
  labs(x = NULL, y = "Estimated fraction (NNLS)", fill = "cell type",
       title = "A  Bulk deconvolution using the D6 single-cell signature (GSE241715, n=21)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 8),
        plot.title = element_text(face = "bold"))

ord <- grp[group == "Bacterial_fungal"][order(-mean_w)]$celltype
grp[, celltype := factor(celltype, levels = rev(ord))]
p2 <- ggplot(grp, aes(x = celltype, y = mean_w, fill = group)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  coord_flip() +
  labs(x = NULL, y = "Mean estimated fraction", fill = NULL,
       title = "B  Cell-type composition across bulk groups") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")

top_cor <- cors[abs(rho) > 0.45][order(gene, -abs(rho))][, head(.SD, 2), by = gene]
p3 <- ggplot(top_cor, aes(x = reorder(paste(gene, celltype, sep = " ~ "), rho), y = rho, fill = rho > 0)) +
  geom_col(width = 0.7) + coord_flip() +
  scale_fill_manual(values = c(`TRUE` = "#C44E52", `FALSE` = "#4C72B0"), guide = "none") +
  geom_hline(yintercept = 0, linewidth = 0.35) +
  labs(x = NULL, y = "Spearman rho (gene expression vs estimated fraction)",
       title = "C  Which cell type explains core-gene expression?") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"))

suppressPackageStartupMessages(library(patchwork))
p <- p1 / (p2 | p3) + plot_layout(heights = c(1, 1.1))
ggsave(file.path(fig, "D7_Fig1_deconvolution.png"), p, width = 15, height = 10, dpi = 300)
ggsave(file.path(fig, "D7_Fig1_deconvolution.pdf"), p, width = 15, height = 10)

cat("\nD7-01 done\n")