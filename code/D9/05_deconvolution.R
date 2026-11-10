# ============================================================
# D9-05  待办 10.3-(3)：正向解卷积 —— 把"上皮效应"落到细胞类型上
#
#   问题：hkManifoldAlignment 给的是"基因层面扰动谱"，不告诉你扰动发生在
#   哪种细胞里。总部报告 Fig9/Fig10 的"epithelial arm"依赖手工定义 18 基因。
#   本脚本用 D6 单细胞的细胞类型特异表达谱做**非负最小二乘（NNLS）解卷积**，
#   把每个 KO 的上皮扰动向量投射到 11 个细胞类型的"表达签名"上，
#   得到"这个敲除主要扰动哪一类细胞的转录组"。
#
#   方法：
#     S  : c x g 细胞类型签名（c=11 细胞类型，g=947 基因，logCPM 均值，行内 Z）
#     y_i: g 维扰动向量（第 i 个 KO 的 signed kernel effect）
#     求解 min ||S' w - y_i|| s.t. w >= 0  →  w 为该 KO 对各细胞类型的扰动权重
#   校验：用 D6 的 AUCell / celltype 富集做外部一致性检查。
#
#   产出：
#     results/D9_deconv_signature.tsv
#     results/D9_deconv_ko_celltype_weights.tsv
#     results/D9_deconv_summary.tsv
#     figures/main/D9_Fig8_deconv.{png,pdf}
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(Seurat); library(Matrix)
  library(quadprog); library(ggplot2); library(scTenifoldKnk)
})
OUT <- "D:/A生信分析/虚拟敲除/第9天执行结果"
D8  <- "D:/A生信分析/虚拟敲除/第8天执行结果"
RES <- file.path(OUT, "results"); FIG <- file.path(OUT, "figures", "main")

CORE <- c("Serpine1","Tmsb4x","P4hb","S100a13","Sod1","Prdx2","Vdac1","Capg","Gja1","Xiap","Vim")

## ---------- 1. 细胞类型签名 ----------
WT  <- as.matrix(readRDS(file.path(D8, "data", "D8_WT_network.rds")))
diag(WT) <- 0
genes <- rownames(WT)

sc <- readRDS("D:/A生信分析/虚拟敲除/第6天执行结果/data/D6_seurat_GSE275678_final.rds")
sc <- JoinLayers(sc)
cnt <- GetAssayData(sc, layer = "counts")[genes, ]
ct  <- as.character(sc$celltype)
cat("细胞类型:\n"); print(table(ct))

# 每个细胞类型的 logCPM 均值
xc <- as(cnt, "dgCMatrix")
lib <- Matrix::colSums(xc)
lcp <- log1p(t(t(xc) / lib * 1e6))
sig <- sapply(split(seq_along(ct), ct), function(idx) Matrix::rowMeans(lcp[, idx, drop = FALSE]))
colnames(sig) <- names(split(seq_along(ct), ct))
rownames(sig) <- genes
cat("\n签名矩阵:", dim(sig), "(基因 x 细胞类型)\n")
fwrite(as.data.table(sig, keep.rownames = "gene"), file.path(RES, "D9_deconv_signature.tsv"), sep = "\t")
rm(sc, cnt, xc, lcp); gc()

## ---------- 2. NNLS：把每个 KO 的扰动向量投到细胞类型上 ----------
H <- scTenifoldKnk::heatKernel(WT, t = 10)
qcX <- readRDS(file.path(D8, "data", "D8_knk_qcCounts.rds"))
S <- scTenifoldKnk:::.logCPM(as(qcX, "dgCMatrix"), genes)
kv <- function(gk) as.numeric(crossprod(-S$mu[gk]/S$sd[gk], H[gk,,drop=FALSE])) * S$sd

# 签名做列内 Z（让各细胞类型可比）：sig 是 947 x 11（基因 × 细胞类型）
# Sz 需转成 11 x 947，使行=细胞类型、列=基因，才能与非负权重对齐
Sz <- scale(sig)                            # 947 x 11：每列（细胞类型）Z
Sz <- t(Sz)                                 # 11 x 947
Sz[is.na(Sz)] <- 0
cat("Sz:", dim(Sz), "(细胞类型 x 基因)\n")
Dmat <- 1e-6 * diag(nrow(Sz))              # 轻微岭正则，稳定
solve_nnls <- function(y) {
  # 单向量求解：min ||G w - y||, w >= 0，其中 G = t(Sz)（947 x 11）
  # 法方程： (Sz Sz^T + lambda I) w = Sz y     —— 11 x 11 系统
  y <- as.numeric(y)
  d <- tcrossprod(Sz) + Dmat               # 11 x 11
  b <- as.numeric(Sz %*% y)                # 11
  tryCatch(solve.QP(d, b, Amat = diag(11), bvec = rep(0, 11), meq = 0)$solution,
           error = function(e) rep(NA_real_, 11))
}

ko_list <- genes
# 一次算出全部 947 个 KO 的扰动向量（g x 947），再逐个 NNLS
Ymat <- vapply(ko_list, function(g) as.numeric(kv(g)), numeric(length(genes)))  # 947 x 947
cat("Ymat:", dim(Ymat), "\n")
W <- t(vapply(seq_len(ncol(Ymat)), function(j) solve_nnls(Ymat[, j]),
              numeric(nrow(Sz))))          # 947 x 11
cat("W:", dim(W), " Sz rownames:", paste(rownames(Sz), collapse=","), "\n")
dimnames(W) <- list(ko_list, colnames(sig))
W[is.na(W)] <- 0
Wn <- W / rowSums(W)                       # 归一为比例
fwrite(as.data.table(Wn, keep.rownames = "ko_gene"),
       file.path(RES, "D9_deconv_ko_celltype_weights.tsv"), sep = "\t")

## ---------- 3. 汇总 ----------
summ <- data.table(celltype = colnames(Wn),
                   mean_weight = colMeans(Wn),
                   core_mean = colMeans(Wn[intersect(CORE, rownames(Wn)), , drop = FALSE]),
                   tmsb4x = Wn["Tmsb4x", ])
setorder(summ, -mean_weight)
fwrite(summ, file.path(RES, "D9_deconv_summary.tsv"), sep = "\t")
cat("\n=== 平均解卷积权重（947 KO）===\n"); print(summ)
cat("\n=== 核心基因在细胞类型上的权重 ===\n")
print(round(Wn[intersect(CORE, rownames(Wn)),
               c("Epithelium_Basal","Epithelium_Wing","Epithelium_Superficial",
                 "Epithelium_Proliferating","Endothelium","Endothelium_Lymphatic",
                 "Macrophage","Neutrophil","Fibroblast")], 3))

## ---------- 4. 图 D9-Fig8：核心基因 × 细胞类型的权重热图 ----------
pl <- melt(as.data.table(Wn[intersect(CORE, rownames(Wn)), , drop = FALSE],
                         keep.rownames = "ko_gene"),
           id.vars = "ko_gene", variable.name = "celltype", value.name = "weight")
pl[, ko_gene := factor(ko_gene, levels = intersect(CORE, rownames(Wn)))]
pl[, celltype := factor(celltype, levels = rev(summ$celltype))]
pl[, weight := as.numeric(weight)]
pl[, fillval := ifelse(weight < 0.005, NA_real_, weight)]
p8 <- ggplot(pl, aes(x = ko_gene, y = celltype)) +
  geom_tile(aes(fill = fillval), colour = "grey90", linewidth = 0.3) +
  geom_point(data = pl[weight < 0.005], aes(x = ko_gene, y = celltype),
             shape = 4, size = 0.9, colour = "grey70", stroke = 0.35) +
  geom_text(aes(label = ifelse(weight >= 0.05, sprintf("%.2f", weight), "")),
            size = 1.9, colour = "grey15") +
  scale_fill_gradientn(colours = c("#FFFFFF", "#FDDBC7", "#F4A582", "#D6604D", "#8B1A1A"),
                       limits = c(0, 1), na.value = "#FAFAFA",
                       name = "deconvolution\nweight") +
  labs(x = NULL, y = NULL,
       title = "Forward deconvolution: which cell type does each knockout perturb?",
       subtitle = "Non-negative least squares of each perturbation vector onto 11 cell-type signatures (947-gene network); weights row-normalised to 1") +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.title = element_text(size = 9.5),
        plot.subtitle = element_text(size = 6.6, colour = "grey30"))
ggsave(file.path(FIG, "D9_Fig8_deconv.png"), p8, width = 7.4, height = 4.0, dpi = 300, bg = "white")
ggsave(file.path(FIG, "D9_Fig8_deconv.pdf"), p8, width = 7.4, height = 4.0, bg = "white")
cat("\n[saved] D9_Fig8_deconv\n")
cat("\n=== DONE 05 ===\n")