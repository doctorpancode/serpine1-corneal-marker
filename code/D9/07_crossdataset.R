# ============================================================
# D9-07  待办 10.3-(4)：GSE326016 / GSE319836 再分析
#
#   目的：把"上皮程序被扰动"的结论，从 GSE275678 的虚拟敲除
#        推到另外两个独立的单细胞/单核图谱上做交叉验证。
#
#   - GSE326016：糖尿病角膜缘图谱（人）
#   - GSE319836：人眼表双潜能基底细胞（人）
#
#   三个层面：
#     (a) 若 GEO 上有可在本地处理的矩阵（10x MTX 或已处理 tsv），直接下载重算；
#     (b) 若没有，用 D1 已经解析过的 pseudobulk 代理（proxy_pseudobulk）做基因级比较；
#     (c) 输出"三数据集 × 10 核心基因"的方向一致性矩阵，作为 Fig10。
#
#   诚实性说明写在输出表 note 列：哪一层数据是真的单细胞、哪一层是代理。
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(Matrix)
})
OUT <- "D:/A生信分析/虚拟敲除/第9天执行结果"
RES <- file.path(OUT, "results"); FIG <- file.path(OUT, "figures", "main")
DAT <- file.path(OUT, "data")

core <- c("Serpine1", "Tmsb4x", "P4hb", "S100a13", "Sod1", "Prdx2",
          "Vdac1", "Capg", "Gja1", "Xiap")
epi <- c("Krt5", "Krt12", "Krt14", "Krt15", "Krt17", "Perp", "Aqp5",
         "Gsta4", "Upk3bl", "Upk1b", "Calml3", "Sfn", "Krt6a", "Krt6b",
         "Lypd2", "Tmprss11d", "S100a11", "Anxa2")

## ---------- 1. 盘点 D1 已解析的 GSE326016 / GSE319836 ----------
idx_path <- "D:/A生信分析/虚拟敲除/D1执行结果/第一天执行结果/data/processed/D1_dataset_index.tsv"
idx <- fread(idx_path)
cat("D1 数据集索引:\n"); print(idx)
targets <- c("GSE326016", "GSE319836")
sub <- idx[grepl(paste(targets, collapse = "|"), dataset)]
print(sub)

proc_dir <- "D:/A生信分析/虚拟敲除/D1执行结果/第一天执行结果/data/processed"
files <- list.files(proc_dir, full.names = TRUE)
tf <- files[grepl("326016|319836", files)]
cat("\n本地已有文件:\n"); print(basename(tf))

## ---------- 2. 逐数据集读取并算核心基因 + 上皮程序 ----------
read_ds <- function(path) {
  if (grepl("\\.tsv$", path)) fread(path)
  else if (grepl("\\.rds$", path)) as.data.table(readRDS(path))
  else NULL
}

rows <- list()
for (p in tf) {
  d <- tryCatch(read_ds(p), error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) next
  cat("\n--", basename(p), "dim:", paste(dim(d), collapse = "x"), "cols:", paste(head(colnames(d), 6), collapse = ","), "\n")
  rows[[basename(p)]] <- d
}

## ---------- 3. 若无单细胞矩阵，用 D6 的 pseudobulk 代理做三数据集方向矩阵 ----------
# D6 已产出：D6_cross_dataset_direction.tsv（GSE275678 + GSE294213 的 10 基因方向）
cd <- fread("D:/A生信分析/虚拟敲除/第6天执行结果/results/D6_cross_dataset_direction.tsv")
cat("\nD6 已有的跨数据集方向表:\n"); print(cd)

dir_w <- dcast(cd, gene ~ dataset, value.var = "direction")
setcolorder(dir_w, c("gene", "GSE275678", "GSE294213"))
fwrite(dir_w, file.path(RES, "D9_crossdataset_direction_matrix.tsv"), sep = "\t")

## ---------- 4. 上皮程序在各数据集的可用性核查 ----------
# 检查 18 个上皮基因在 GSE294213（人）里是否可测
g294 <- fread("D:/A生信分析/虚拟敲除/第6天执行结果/results/D6_GSE294213_top_markers.tsv")
cat("\nGSE294213 top markers 前 10:\n"); print(head(g294, 10))

## ---------- 5. 三数据集一致性汇总 ----------
# 已知：GSE275678（D6 单细胞 + D7/D8 虚拟敲除）、GSE294213（D6 跨数据集）
# 对 GSE326016 / GSE319836 用 GEO 元数据 + 本研究结论做"文献层面"的交叉对照，
# 明确标注为 literature-level，不与实测混用。
lit <- data.table(
  dataset = c("GSE275678", "GSE294213", "GSE326016", "GSE319836"),
  species = c("Mouse", "Human", "Human", "Human"),
  assay = c("scRNA-seq", "scRNA-seq", "scRNA-seq / snRNA-seq", "snRNA-seq"),
  condition = c("Alkali burn vs uninjured", "Suture-induced CNV vs normal",
                "Diabetic limbal landscape", "Human ocular surface bipotent basal cells"),
  evidence_level = c("measured_in_project", "measured_in_project",
                     "literature_only", "literature_only"),
  role_in_study = c("虚拟敲除主数据集", "跨物种验证", "角膜缘上皮身份的外部对照",
                    "基底细胞命运的外部对照"),
  supports_epithelial_program = c("yes (Krt12/Aqp5/Gsta4/Upk3bl perturbed)",
                                  "yes (Capg/Gja1/Prdx2 down, P4hb/S100a13/Serpine1 up)",
                                  "context only", "context only")
)
fwrite(lit, file.path(RES, "D9_crossdataset_evidence_level.tsv"), sep = "\t")
cat("\n=== 证据层级 ===\n"); print(lit)

## ---------- 6. 图：10 核心基因在两个实测数据集里的方向一致性 ----------
# cd 的列是 dataset/gene/Control/Lesion/delta/direction，没有 logFC
w <- dcast(cd, gene ~ dataset, value.var = "delta")
setnames(w, c("gene", "GSE275678", "GSE294213"))
w[, concordant := sign(GSE275678) == sign(GSE294213)]
w[, gene := factor(gene, levels = w[order(GSE275678), gene])]
fwrite(w, file.path(RES, "D9_core_gene_direction_concordance.tsv"), sep = "\t")

pm <- melt(w, id.vars = c("gene", "concordant"), variable.name = "dataset", value.name = "delta")
pm[, dataset := factor(dataset, levels = c("GSE275678", "GSE294213"))]

p <- ggplot(pm, aes(x = gene, y = delta, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.72, linewidth = 0) +
  geom_hline(yintercept = 0, linewidth = 0.35) +
  geom_point(data = pm[concordant == TRUE], aes(x = gene, y = delta),
             inherit.aes = FALSE, shape = 21, size = 0.9, colour = "white",
             fill = "white", alpha = 0.85) +
  scale_fill_manual(values = c(GSE275678 = "#B2182B", GSE294213 = "#4393C3")) +
  labs(x = NULL, y = "expression delta (lesion - control)",
       fill = "dataset",
       title = "Core-gene direction concordance across two lesion datasets",
       subtitle = sprintf("%d / %d genes share the same direction (white dots mark concordant genes)",
                          sum(w$concordant), nrow(w))) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = 6.6, colour = "grey30"))
ggsave(file.path(FIG, "D9_Fig10_crossdataset_concordance.png"), p,
       width = 5.4, height = 3.2, dpi = 300, bg = "white")
ggsave(file.path(FIG, "D9_Fig10_crossdataset_concordance.pdf"), p, width = 5.4, height = 3.2)

cat("\n=== 方向一致率:", sum(w$concordant), "/", nrow(w), "===\n")

## ---------- 7. 真实使用 GSE326016（糖尿病 vs 正常角膜缘）做第三个实测数据集 ----------
# 注意：.rds 读回来的 dgCMatrix **rownames 丢失**（D1 保存时未保留 dimnames），
#       所以必须读同名 .tsv（第一列 SYMBOL，其余为样本 count）。
#   DM19_56 / DM20_01 / DM20_20 = 糖尿病组; N19_02 / N19_52 / N19_55 = 正常组
pb_path <- "D:/A生信分析/虚拟敲除/D1执行结果/第一天执行结果/data/processed/GSE326016_pseudobulk.tsv"
pb <- fread(pb_path)
cn <- setdiff(colnames(pb), "SYMBOL")
mat <- as.matrix(pb[, ..cn]); rownames(mat) <- pb[[1]]
grp <- ifelse(grepl("\\|DM", cn), "DM", "Normal")
cat("\nGSE326016 分组:", paste(cn, grp, sep = "=", collapse = ", "), "\n")
cat("矩阵:", paste(dim(mat), collapse = " x "), "  首行基因:", rownames(mat)[1], "\n")
cat("SERPINE1 可测:", "SERPINE1" %in% rownames(mat), " KRT5 可测:", "KRT5" %in% rownames(mat), "\n")

# CPM -> log2
cpm <- sweep(mat, 2, colSums(mat), "/") * 1e6
l2 <- log2(cpm + 1)
rownames(l2) <- toupper(rownames(l2))   # 人基因名统一大写

# 人-鼠同源映射（只处理本研究涉及的基因）
# 小鼠基因名是 Serpine1 / Tmsb4x 这种 "首字母大写其余小写" 形式，
# 人基因名是全大写 SERPINE1 / TMSB4X。整串 toupper 即可完成映射。
mu2hu_fun <- function(g) toupper(g)

probe <- data.table(mouse = c(core, epi))
probe[, human := vapply(mouse, mu2hu_fun, character(1))]
probe[, present := human %in% rownames(l2)]
cat("可测核心基因:", sum(probe[1:length(core), present]), "/", length(core),
    "；可测上皮程序基因:", sum(probe[(length(core)+1):.N, present]), "/", length(epi), "\n")

test_g <- function(gs, tag) {
  gs <- gs[gs %in% rownames(l2)]
  if (length(gs) == 0) return(data.table(gene_human = character(), DM = numeric(),
                                         Normal = numeric(), delta = numeric(),
                                         t = numeric(), p = numeric(), set = character()))
  rbindlist(lapply(gs, function(g) {
    v <- l2[g, ]
    a <- v[grp == "DM"]; b <- v[grp == "Normal"]
    tt <- tryCatch(t.test(a, b), error = function(e) NULL)
    data.table(gene_human = g, DM = mean(a), Normal = mean(b),
               delta = mean(a) - mean(b),
               t = if (is.null(tt)) NA_real_ else unname(tt$statistic),
               p = if (is.null(tt)) NA_real_ else unname(tt$p.value), set = tag)
  }))
}
r326_core <- test_g(probe[1:length(core), human], "core")
r326_epi  <- test_g(probe[(length(core)+1):.N, human], "epithelial")
r326 <- rbind(r326_core, r326_epi)
r326[, p_adj := p.adjust(p, "BH")]
fwrite(r326, file.path(RES, "D9_GSE326016_core_and_epi.tsv"), sep = "\t")
cat("\n=== GSE326016（糖尿病 vs 正常）核心基因 ===\n"); print(r326_core[order(p)])
cat("\n=== GSE326016 上皮程序基因（按 p 排序前 10）===\n"); print(head(r326_epi[order(p)], 10))
cat("\n共识方向：核心基因里 DM 上调/下调数 =",
    nrow(r326_core[delta > 0]), "/", nrow(r326_core[delta < 0]), "\n")

## ---------- 8. 上皮程序打分在 GSE326016 里是否整体偏移 ----------
epi_h <- probe[(length(core)+1):.N][present == TRUE, human]
if (length(epi_h) >= 3) {
  z <- t(scale(t(l2[epi_h, ])))
  score <- colMeans(z, na.rm = TRUE)
  sc <- data.table(sample = cn, group = grp, epi_score = as.numeric(score))
  fwrite(sc, file.path(RES, "D9_GSE326016_epi_score.tsv"), sep = "\t")
  cat("\n=== GSE326016 上皮程序打分 ===\n"); print(sc)
  sc_tt <- t.test(epi_score ~ group, data = sc)
  cat("上皮程序打分 t 检验 p =", signif(sc_tt$p.value, 4), "\n")
} else {
  cat("\n上皮程序基因可测数不足（", length(epi_h), "），跳过打分。\n")
}

## ---------- 9. 汇总：三数据集方向矩阵 ----------
g326 <- merge(r326_core[, .(gene_human, GSE326016 = sign(delta))],
              unique(probe[, .(mouse, human)]), by.x = "gene_human", by.y = "human")
w3 <- merge(w, g326[, .(gene = mouse, GSE326016)], by = "gene", all.x = TRUE)
w3[, concordant_all := (sign(GSE275678) == sign(GSE294213)) &
                      (is.na(GSE326016) | sign(GSE275678) == GSE326016)]
fwrite(w3, file.path(RES, "D9_core_gene_direction_concordance_3ds.tsv"), sep = "\t")
cat("\n=== 三数据集方向矩阵 ===\n"); print(w3)

if (nrow(r326_core) > 0) {
  pm2 <- rbind(
    data.table(gene = probe[1:length(core)][human %in% r326_core$gene_human][match(r326_core$gene_human, human), mouse],
               dataset = "GSE326016", delta = r326_core$delta, p = r326_core$p)
  )
  p2 <- ggplot(pm2, aes(x = gene, y = delta)) +
    geom_col(fill = "#7B3294", width = 0.68, linewidth = 0) +
    geom_hline(yintercept = 0, linewidth = 0.35) +
    labs(x = NULL, y = expression(log[2]~"CPM delta (diabetic - normal)"),
         title = "GSE326016: core genes in diabetic vs normal limbus") +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          panel.grid.minor = element_blank())
  ggsave(file.path(FIG, "D9_Fig10b_GSE326016_core.png"), p2, width = 4.4, height = 3.0, dpi = 300, bg = "white")
  ggsave(file.path(FIG, "D9_Fig10b_GSE326016_core.pdf"), p2, width = 4.4, height = 3.0)
}

cat("\n=== DONE 07 ===\n")