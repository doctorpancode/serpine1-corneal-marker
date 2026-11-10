# ============================================================
# D9-01  待办 10.1-(1)：补一个"炎症性人角膜全层"验证队列  (v2)
#
#   问题：总报告唯一的人验证集是 GSE204791（圆锥角膜 KC），
#         KC 是非炎症性变薄，SERPINE1 在 KC 基质里**下调**，
#         与炎症性 CNV 的**上调**方向冲突 → 审稿人必打。
#
#   方案：补 GSE58291（Illumina HumanHT-12 v4，GPL10558，30 例）
#         Control 12 / 真菌性角膜炎 FK 8 / 细菌性角膜炎 BK 7
#         + 作者已剔除 3 个离群样本（M50/M73/M80）。
#
#   关键工程决策（踩坑记录）：
#     1) 原始文件 GSE58291_Non-normalized_data_with_Gene_ID.txt.gz 的
#        header 与数据行**列数不等**（header 270 列、数据行 286 列，
#        中间少了 16 个字段），用列名无法对齐；
#        但**前 2 + 30*9 = 272 列的位置是对齐的**，故按位置取
#        `2 + 9*30` 的索引，用 header 解析出的样本名。
#     2) 作者已在 series 级提供三份 VST+quantile 归一化矩阵
#        （Bacterial_v_Control / Fungal_v_Control / Bacterial_v_Fungal），
#        样本列名是 **M23/T4/M6…** 这类内部编号。**这才是可用的分析对象。**
#     3) 内部编号 → GSM 的映射官方未提供，故本脚本**不猜映射**，
#        只用样本名做"跨比较的内在对照"：追踪同一个 M 编号在三份
#        子矩阵中的表达，即可判定它属于哪一组（分组内应稳定、跨组应
#        大幅位移）。
#
#   产出：
#     results/D9_GSE58291_grouping_recovery.tsv   分组恢复的证据
#     results/D9_GSE58291_sample_group.tsv        最终样本分组
#     results/D9_core_genes_GSE58291.tsv          10 核心基因三组表达 + 统计
#     results/D9_GSE58291_limma_*.tsv             差异分析
#     results/D9_crossplatform_consistency.tsv    与 GSE241715 方向一致性
#     figures/main/D9_Fig1*, D9_Fig2*
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(limma); library(ggplot2)
})
OUT  <- "D:/A生信分析/虚拟敲除/第9天执行结果"
DATA <- file.path(OUT, "data"); RES <- file.path(OUT, "results")
FIG  <- file.path(OUT, "figures", "main")
for (d in c(RES, FIG)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

CORE <- c("SERPINE1","TMSB4X","P4HB","S100A13","SOD1","PRDX2","VDAC1","CAPG","GJA1","XIAP")

## ------------------------------------------------------------------
## 1. 探针 ID → SYMBOL（来自原始文件的按位置读取）
## ------------------------------------------------------------------
fraw <- file.path(DATA, "GSE58291_expr.tsv")
con <- file(fraw, "r"); hdr_line <- readLines(con, n = 1); close(con)
hdr <- strsplit(hdr_line, "\t", fixed = TRUE)[[1]]
n_avg <- sum(grepl("\\.AVG_Signal$", hdr))
stopifnot(n_avg == 30)
avg_idx  <- which(grepl("\\.AVG_Signal$", hdr))     # 30 个位置
samp_name <- sub("^[MT]\\s*", "", sub("\\.AVG_Signal$", "", hdr[avg_idx]))
cat("AVG_Signal 样本名（内部编号）:", paste(samp_name, collapse = ", "), "\n")

raw <- fread(fraw, sep = "\t", header = FALSE, skip = 1, select = c(1, 2, avg_idx),
             col.names = c("probe", "symbol", samp_name), showProgress = FALSE)
cat("原始行数:", nrow(raw), "\n")
probe2sym <- unique(raw[, .(probe, symbol = toupper(trimws(symbol)))])

## ------------------------------------------------------------------
## 2. 读三份作者归一化矩阵（VST + quantile，log2 尺度）
## ------------------------------------------------------------------
read_norm <- function(f) {
  d <- fread(file.path(DATA, f), sep = "\t", header = TRUE, showProgress = FALSE)
  setnames(d, 1, "probe")
  keep <- !grepl("\\.Detection\\.Pval$", names(d)) & names(d) != "probe"
  d <- d[, c("probe", names(d)[keep]), with = FALSE]
  m <- as.matrix(d[, -1]); storage.mode(m) <- "double"
  rownames(m) <- d$probe
  m
}
m_bc <- read_norm("GSE58291_BK_vs_Control_normalized.tsv")
m_fc <- read_norm("GSE58291_FK_vs_Control_normalized.tsv")
m_bf <- read_norm("GSE58291_BK_vs_FK_normalized.tsv")
cat("BK_vs_Ctrl:", dim(m_bc), " cols:", paste(colnames(m_bc), collapse = ","), "\n")
cat("FK_vs_Ctrl:", dim(m_fc), " cols:", paste(colnames(m_fc), collapse = ","), "\n")
cat("BK_vs_FK  :", dim(m_bf), " cols:", paste(colnames(m_bf), collapse = ","), "\n")

## ------------------------------------------------------------------
## 3. 分组恢复
##    Control 样本 = M12,M23,M26,M30,M31,M33,M34,M43,M44,M45,M49,M6,M61,M66
##    （在 BK_vs_Ctrl 与 FK_vs_Ctrl 都出现的那批）
##    BK = 只在 BK_vs_Ctrl 出现；FK = 只在 FK_vs_Ctrl 出现
## ------------------------------------------------------------------
in_bc <- colnames(m_bc); in_fc <- colnames(m_fc)
ctrl <- intersect(in_bc, in_fc)
bk   <- setdiff(in_bc, ctrl)
fk   <- setdiff(in_fc, ctrl)
cat("\n交集(Control):", length(ctrl), paste(sort(ctrl), collapse = ","), "\n")
cat("仅 BK 组:", length(bk), paste(sort(bk), collapse = ","), "\n")
cat("仅 FK 组:", length(fk), paste(sort(fk), collapse = ","), "\n")

grp_tab <- data.table(sample = c(ctrl, fk, bk),
                      group  = c(rep("Control", length(ctrl)),
                                 rep("FK", length(fk)),
                                 rep("BK", length(bk))))
# 注：样本名来自归一化矩阵（M12/T4…），原始文件的 AVG_Signal 列名是 "M 12"/"T 4"
# （带空格），两者必须 stripspace 后比对；此处只做归一化矩阵内部的完整性检查。
norm_samples <- c(in_bc, in_fc)
stopifnot(all(c(ctrl, fk, bk) %in% norm_samples))
fwrite(grp_tab, file.path(RES, "D9_GSE58291_sample_group.tsv"), sep = "\t")

## 3b. 用"同一 M 编号在三份子矩阵中的表达"验证分组（不依赖任何外部映射）
common_probes <- Reduce(intersect, list(rownames(m_bc), rownames(m_fc), rownames(m_bf)))
cat("三份矩阵共有探针:", length(common_probes), "\n")
in_bf <- colnames(m_bf)
evid <- rbindlist(lapply(c(ctrl, fk, bk), function(s) {
  data.table(sample = s, group = grp_tab[sample == s, group],
             in_BKvCtrl = s %in% in_bc, in_FKvCtrl = s %in% in_fc, in_BKvFK = s %in% in_bf,
             mean_expr_BKvCtrl = if (s %in% in_bc) mean(m_bc[common_probes, s]) else NA_real_,
             mean_expr_FKvCtrl = if (s %in% in_fc) mean(m_fc[common_probes, s]) else NA_real_)
}))
fwrite(evid, file.path(RES, "D9_GSE58291_grouping_recovery.tsv"), sep = "\t")

## ------------------------------------------------------------------
## 4. 合并成一张矩阵（取三份的并集探针，缺失用 NA）
## ------------------------------------------------------------------
allsamp <- c(ctrl, fk, bk)
M <- matrix(NA_real_, nrow = length(common_probes), ncol = length(allsamp),
            dimnames = list(common_probes, allsamp))
M[, intersect(ctrl, in_bc)] <- m_bc[common_probes, intersect(ctrl, in_bc)]
missing_ctrl <- setdiff(ctrl, in_bc)
if (length(missing_ctrl)) M[, missing_ctrl] <- m_fc[common_probes, missing_ctrl]
M[, bk] <- m_bc[common_probes, bk]
M[, fk] <- m_fc[common_probes, fk]
cat("合并矩阵:", dim(M), " 缺失值:", sum(is.na(M)), "\n")
M <- M[rowSums(is.na(M)) == 0, , drop = FALSE]
cat("完整观测探针:", nrow(M), "\n")

## 探针 → 基因（同基因多探针取方差最大的探针）
p2s <- probe2sym[probe %in% rownames(M)]
p2s <- p2s[!is.na(symbol) & symbol != "" & symbol != "NA"]
v <- apply(M[p2s$probe, , drop = FALSE], 1, var)
ord <- order(v, decreasing = TRUE)
p2s <- p2s[ord][!duplicated(p2s$symbol[ord])]
Mg <- M[p2s$probe, , drop = FALSE]
rownames(Mg) <- p2s$symbol
cat("基因级矩阵:", dim(Mg), "\n")

grp <- factor(grp_tab[match(colnames(Mg), sample), group], levels = c("Control","FK","BK"))
stopifnot(!any(is.na(grp)))
cat("分组:", paste(names(table(grp)), table(grp), sep = "=", collapse = ", "), "\n")

## ------------------------------------------------------------------
## 5. limma（BK vs Control, FK vs Control）
## ------------------------------------------------------------------
design <- model.matrix(~ 0 + grp); colnames(design) <- levels(grp)
fit  <- lmFit(Mg, design)
cont <- makeContrasts(BK_vs_Ctrl = BK - Control, FK_vs_Ctrl = FK - Control,
                      BK_vs_FK = BK - FK, levels = design)
fit2 <- eBayes(contrasts.fit(fit, cont))
save_limma <- function(coef, name) {
  r <- as.data.table(topTable(fit2, coef = coef, number = Inf, sort.by = "P"), keep.rownames = "gene")
  setnames(r, c("gene","logFC","AveExpr","t","P.Value","adj.P.Val","B"))
  fwrite(r, file.path(RES, paste0("D9_GSE58291_limma_", name, ".tsv")), sep = "\t")
  r
}
r_bk <- save_limma("BK_vs_Ctrl", "BK_vs_Control")
r_fk <- save_limma("FK_vs_Ctrl", "FK_vs_Control")
r_bf <- save_limma("BK_vs_FK",   "BK_vs_FK")
cat("\nBK vs Control: DEG(FDR<0.05) =", r_bk[adj.P.Val < 0.05, .N], "\n")
cat("FK vs Control: DEG(FDR<0.05) =", r_fk[adj.P.Val < 0.05, .N], "\n")

## ------------------------------------------------------------------
## 6. 10 个核心基因
## ------------------------------------------------------------------
long <- melt(as.data.table(Mg[intersect(CORE, rownames(Mg)), , drop = FALSE], keep.rownames = "gene"),
             id.vars = "gene", variable.name = "sample", value.name = "expr")
long[, group := grp[match(sample, colnames(Mg))]]
long[, group := factor(group, levels = c("Control","FK","BK"))]

stat <- long[, .(n = .N, mean = mean(expr), median = median(expr), sd = sd(expr)), by = .(gene, group)]
wide <- dcast(stat, gene ~ group, value.var = c("mean","median","n"))
wide <- merge(wide, r_bk[, .(gene, BK_logFC = logFC, BK_p = P.Value, BK_adjP = adj.P.Val)], by = "gene", all.x = TRUE)
wide <- merge(wide, r_fk[, .(gene, FK_logFC = logFC, FK_p = P.Value, FK_adjP = adj.P.Val)], by = "gene", all.x = TRUE)
kw <- long[, {
  kt <- kruskal.test(expr ~ group)
  ct <- suppressWarnings(cor.test(expr, as.integer(group), method = "spearman", exact = FALSE))
  .(kw_p = kt$p.value, trend_rho = unname(ct$estimate), trend_p = ct$p.value)
}, by = gene]
wide <- merge(wide, kw, by = "gene")
setcolorder(wide, c("gene","mean_Control","mean_FK","mean_BK",
                    "BK_logFC","BK_p","BK_adjP","FK_logFC","FK_p","FK_adjP",
                    "kw_p","trend_rho","trend_p"))
setorder(wide, -mean_BK)
fwrite(wide, file.path(RES, "D9_core_genes_GSE58291.tsv"), sep = "\t")

cat("\n=== 10 核心基因 @ GSE58291（VST log2）===\n")
print(wide[, .(gene, mean_Control = round(mean_Control,2), mean_FK = round(mean_FK,2),
               mean_BK = round(mean_BK,2), BK_logFC = round(BK_logFC,2),
               BK_adjP = signif(BK_adjP,3), trend_rho = round(trend_rho,2))])

## ------------------------------------------------------------------
## 7. 与 GSE241715 的跨平台方向一致性
## ------------------------------------------------------------------
d3 <- "D:/A生信分析/虚拟敲除/第3天执行结果/results/D3_limma_deg.tsv"
if (file.exists(d3)) {
  deg <- fread(d3)
  cat("\nD3 limma 列名:", paste(names(deg), collapse = ", "), "\n")
  if (!"gene" %in% names(deg)) setnames(deg, 1, "gene")
  anchor <- intersect(CORE, deg$gene)
  cmp <- merge(wide[gene %in% anchor, .(gene, BK_logFC, FK_logFC, BK_adjP)],
               deg[gene %in% anchor, .(gene, R715_logFC = logFC, R715_adjP = adjP)], by = "gene")
  cmp[, consistent_BK := sign(BK_logFC) == sign(R715_logFC)]
  cmp[, consistent_FK := sign(FK_logFC) == sign(R715_logFC)]
  fwrite(cmp, file.path(RES, "D9_crossplatform_consistency.tsv"), sep = "\t")
  cat("\n=== 跨平台方向一致性（GSE58291 vs GSE241715）===\n"); print(cmp)
  cat("BK 一致:", sum(cmp$consistent_BK), "/", nrow(cmp),
      "| FK 一致:", sum(cmp$consistent_FK), "/", nrow(cmp), "\n")
}

## ------------------------------------------------------------------
## 8. 图 D9-Fig1（小面箱线）与 D9-Fig2（梯度折线）
## ------------------------------------------------------------------
pdt <- copy(long[gene %in% CORE])
pdt <- merge(pdt, unique(wide[, .(gene, lab = ifelse(BK_adjP < 0.001, "***",
                                       ifelse(BK_adjP < 0.01, "**",
                                       ifelse(BK_adjP < 0.05, "*", "ns"))))]), by = "gene", all.x = TRUE)
pdt[, gene := factor(gene, levels = CORE)]

p1 <- ggplot(pdt, aes(x = group, y = expr, fill = group)) +
  geom_boxplot(width = 0.62, outlier.size = 0.25, linewidth = 0.28) +
  geom_jitter(width = 0.1, size = 0.7, alpha = 0.55, colour = "grey20") +
  geom_text(data = unique(pdt[, .(gene, group = factor("BK", levels = levels(group)), lab)]),
            aes(x = group, y = Inf, label = lab), vjust = 1.5, size = 2.6, inherit.aes = FALSE) +
  facet_wrap(~ gene, nrow = 2, scales = "free_y") +
  scale_fill_manual(values = c(Control = "#8DA0CB", FK = "#FC8D62", BK = "#B2182B")) +
  labs(x = NULL, y = expression(log[2]~"VST signal"),
       title = "Human inflamed cornea (GSE58291): core genes vs injury severity",
       subtitle = "Control n=12 (cadaveric) | FK n=8 | BK n=7 (after removing 3 author-flagged outliers)   *** FDR < 0.001 (BK vs Control)") +
  theme_bw(base_size = 9) +
  theme(legend.position = "none", strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text = element_text(face = "bold", size = 8),
        panel.grid.minor = element_blank(), axis.text.x = element_text(size = 7))
ggsave(file.path(FIG, "D9_Fig1_GSE58291_core_genes.png"), p1, width = 9, height = 4.4, dpi = 300)
ggsave(file.path(FIG, "D9_Fig1_GSE58291_core_genes.pdf"), p1, width = 9, height = 4.4)
cat("\n[saved] D9_Fig1\n")

g3 <- melt(wide[, .(gene, Control = mean_Control, FK = mean_FK, BK = mean_BK)],
           id.vars = "gene", variable.name = "group", value.name = "mean")
g3[, group := factor(group, levels = c("Control","FK","BK"))]
g3[, gene := factor(gene, levels = wide[order(mean_BK)]$gene)]
p2 <- ggplot(g3, aes(x = group, y = mean, group = gene, colour = gene)) +
  geom_line(linewidth = 0.55, alpha = 0.85) + geom_point(size = 1.7) +
  scale_colour_manual(values = c(SERPINE1 = "#B2182B", TMSB4X = "#2166AC",
                                 P4HB = "#F4A582", S100A13 = "#92C5DE",
                                 SOD1 = "#66A61E", PRDX2 = "#E6AB02",
                                 VDAC1 = "#7570B3", CAPG = "#A6761D",
                                 GJA1 = "#1B9E77", XIAP = "#666666")) +
  labs(x = "Injury severity", y = expression("mean " * log[2] * " VST signal"),
       title = "Dose-response across injury severity (GSE58291)", colour = NULL) +
  theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank())
ggsave(file.path(FIG, "D9_Fig2_GSE58291_gradient.png"), p2, width = 5.2, height = 3.8, dpi = 300)
ggsave(file.path(FIG, "D9_Fig2_GSE58291_gradient.pdf"), p2, width = 5.2, height = 3.8)
cat("[saved] D9_Fig2\n")

cat("\n=== DONE 01 ===\n")