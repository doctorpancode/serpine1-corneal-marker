# =============================================================
# D2-02  StemnessScore (ssGSEA) + WGCNA（限血管生成/角化/ECM 基因池）
# 输入:
#   D1: D:/A生信分析/虚拟敲除/D1执行结果/D1_CNV_角膜新生血管_数据分析/data/processed/
#        GSE241715_expression.rds (ExpressionSet 16463 x 21, log2(CPM+1))
#        GSE204791_expression.rds (ExpressionSet 58339 x 31, probe x sample)
#        GSE204791_sample_sheet.tsv
#   D2: ../data/stemness_genesets.tsv
# 输出:
#   data/D2_expression_analysis.rds        (合并后的分析对象: QC 后的 gene x sample 矩阵 + 表型)
#   data/stemness_score.tsv               (每样本 ssGSEA 分数)
#   data/wgcna_gene_pool.tsv              (限维后的基因池及来源通路)
#   data/wgcna_module_assignment.tsv
#   data/wgcna_module_trait_corr.tsv
#   results/wgcna_softthreshold.tsv
#   results/wgcna_summary.tsv
#   figures/D2_Fig1_softthreshold.pdf/.png
#   figures/D2_Fig2_module_trait_heatmap.pdf/.png
#   figures/D2_Fig3_stemness_boxplot.pdf/.png
#   figures/D2_Fig4_dendrogram.pdf/.png
# =============================================================
set.seed(20261109)
suppressPackageStartupMessages({
  library(Biobase); library(GSVA); library(WGCNA); library(data.table)
  library(ggplot2); library(patchwork); library(pheatmap); library(msigdbr)
})
options(stringsAsFactors = FALSE)
allowWGCNAThreads(nThreads = 8)

d1  <- "D:/A生信分析/虚拟敲除/D1执行结果/D1_CNV_角膜新生血管_数据分析/data/processed"
out <- "D:/A生信分析/虚拟敲除/第2天执行结果"
dir.create(file.path(out, "data"),    recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)

log_con <- file(file.path(out, "logs", "02_stemness_wgcna.log"), open = "wt")
sink(log_con, split = TRUE)
cat("=== D2-02 start:", as.character(Sys.time()), "===\n")

# -------------------------------------------------------------
# 1. GSE241715（人角膜 bulk RNA-seq，21 例）——主分析队列
# -------------------------------------------------------------
eset <- readRDS(file.path(d1, "GSE241715_expression.rds"))
expr <- exprs(eset)
pd   <- pData(eset)
cat("GSE241715:", nrow(expr), "genes x", ncol(expr), "samples\n")
print(table(pd$group))

# 分组：病变 vs 正常（二分类），保留病因分型
grp <- factor(pd$group)
binary <- factor(ifelse(grp == "Control", "Normal", "Lesion"),
                 levels = c("Normal", "Lesion"))
pd$binary_lesion <- binary
pheno <- data.frame(
  group     = as.character(grp),
  binary    = as.character(binary),
  dataset   = "GSE241715",
  stringsAsFactors = FALSE
)
rownames(pheno) <- colnames(expr)

# -------------------------------------------------------------
# 2. StemnessScore: GSVA ssGSEA
# -------------------------------------------------------------
gs <- fread(file.path(out, "data", "stemness_genesets.tsv"))
gs_list <- list(Corneal_Limbal_Stemness = unique(gs$gene))
gs_list <- lapply(gs_list, function(g) intersect(g, rownames(expr)))
cat("geneset genes present:", length(gs_list[[1]]), "/", length(unique(gs$gene)), "\n")

ssPar <- ssgseaParam(expr, gs_list, normalize = TRUE, minSize = 5, maxSize = 500)
ss <- gsva(ssPar, verbose = FALSE)
stem_score <- as.numeric(ss[1, ])
names(stem_score) <- colnames(ss)

stem_dt <- data.table(sample = names(stem_score),
                      StemnessScore = as.numeric(stem_score),
                      group = pheno[names(stem_score), "group"],
                      binary = pheno[names(stem_score), "binary"])
fwrite(stem_dt, file.path(out, "data", "stemness_score.tsv"), sep = "\t")
cat("\nStemnessScore by group:\n"); print(stem_dt[, .(n = .N, mean = round(mean(StemnessScore), 4),
                                                     sd = round(sd(StemnessScore), 4)), by = group])
wt <- wilcox.test(StemnessScore ~ binary, data = stem_dt)
cat("Wilcoxon Normal vs Lesion p =", signif(wt$p.value, 4), "\n")

fwrite(data.table(test = "Wilcoxon rank-sum Normal vs Lesion",
                  statistic = wt$statistic, p_value = wt$p.value,
                  mean_Normal = mean(stem_dt[binary == "Normal", StemnessScore]),
                  mean_Lesion = mean(stem_dt[binary == "Lesion", StemnessScore])),
       file.path(out, "results", "stemness_group_test.tsv"), sep = "\t")

# -------------------------------------------------------------
# 3. 限定基因池：血管生成 + ECM + 角化（MSigDB + KEGG）
# -------------------------------------------------------------
pool <- character(0)
pool_src <- list()
add_sets <- function(gs_names, collection) {
  db <- msigdbr(species = "Homo sapiens", collection = collection)
  for (nm in gs_names) {
    sub <- db[db$gs_name == nm, ]
    if (nrow(sub) == 0) { cat("  [miss]", nm, "\n"); next }
    pool <<- c(pool, sub$gene_symbol)
    for (g in sub$gene_symbol) pool_src[[g]] <<- paste0(c(pool_src[[g]], nm), collapse = ";")
    cat("  [ok]", nm, nrow(sub), "\n")
  }
}
cat("\n--- building restricted gene pool ---\n")
add_sets(c("GOBP_BLOOD_VESSEL_MORPHOGENESIS",
           "GOBP_BLOOD_VESSEL_REMODELING",
           "GOBP_BLOOD_VESSEL_MATURATION",
           "GOBP_KERATINOCYTE_DIFFERENTIATION",
           "GOBP_KERATINOCYTE_DEVELOPMENT",
           "GOBP_EPIDERMAL_CELL_DIFFERENTIATION",
           "GOBP_EPITHELIUM_DEVELOPMENT",
           "GOBP_CELL_SUBSTRATE_ADHESION",
           "GOBP_RESPONSE_TO_TRANSFORMING_GROWTH_FACTOR_BETA",
           "GOBP_ENDOTHELIAL_CELL_PROLIFERATION",
           "GOBP_ENDOTHELIAL_CELL_MIGRATION"), "C5")
add_sets(c("KEGG_ECM_RECEPTOR_INTERACTION", "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
           "KEGG_TGF_BETA_SIGNALING_PATHWAY", "KEGG_FOCAL_ADHESION",
           "KEGG_VEGF_SIGNALING_PATHWAY", "KEGG_WNT_SIGNALING_PATHWAY"), "C2")
pool <- unique(pool)
pool_present <- intersect(pool, rownames(expr))
cat("pool genes (total / in GSE241715):", length(pool), "/", length(pool_present), "\n")

pool_dt <- data.table(gene = pool_present,
                      sources = sapply(pool_present, function(g) pool_src[[g]]))
# 剔除低表达（<25% 样本未检出水平）与近零方差
keep_expr <- apply(expr[pool_present, ], 1, function(x) sum(x > 0) >= 0.75 * length(x))
pool_dt <- pool_dt[gene %in% names(keep_expr)[keep_expr]]
vars <- apply(expr[pool_dt$gene, ], 1, var)
pool_dt <- pool_dt[gene %in% names(vars)[vars > quantile(vars, 0.05)]]
fwrite(pool_dt, file.path(out, "data", "wgcna_gene_pool.tsv"), sep = "\t")
cat("pool after expression filter:", nrow(pool_dt), "genes\n")

# -------------------------------------------------------------
# 4. WGCNA
# -------------------------------------------------------------
datExpr <- t(expr[pool_dt$gene, ])   # samples x genes
cat("WGCNA input:", nrow(datExpr), "samples x", ncol(datExpr), "genes\n")
gsg <- goodSamplesGenes(datExpr, verbose = 1)
if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
cat("after goodSamplesGenes:", paste(dim(datExpr), collapse = " x "), "\n")

# 样本树 + 表型树（离群样本检查）
sampleTree <- hclust(dist(datExpr), method = "average")
pdf(file.path(out, "figures", "D2_Fig4_dendrogram.pdf"), width = 9, height = 5)
par(cex = 0.7, mar = c(3, 5, 2, 1))
plot(sampleTree, main = "Sample clustering (GSE241715, restricted pool)",
     sub = "", xlab = "", cex.lab = 1.1, cex.axis = 1.0, cex.main = 1.3)
dev.off()
png(file.path(out, "figures", "D2_Fig4_dendrogram.png"), width = 2700, height = 1500, res = 300)
par(cex = 0.7, mar = c(3, 5, 2, 1))
plot(sampleTree, main = "Sample clustering (GSE241715, restricted pool)",
     sub = "", xlab = "", cex.lab = 1.1, cex.axis = 1.0, cex.main = 1.3)
dev.off()

# 软阈值
powers <- c(1:10, 12:20, 22:30)
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed",
                         corFnc = "cor", corOptions = list(use = "p"),
                         verbose = 5)
sft_dt <- data.table(power = sft$powerEstimate, fitInd = NA_real_)
sft_tab <- sft$fitIndices
fwrite(as.data.table(sft_tab), file.path(out, "results", "wgcna_softthreshold.tsv"), sep = "\t")

p1 <- ggplot(sft_tab, aes(Power, -sign(slope) * SFT.R.sq)) +
  geom_text(aes(label = Power), size = 3) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "red") +
  labs(x = "Soft threshold (power)", y = "Scale-free topology fit R^2 (signed)",
       title = "Scale independence") +
  theme_bw(base_size = 12)
p2 <- ggplot(sft_tab, aes(Power, mean.k.)) +
  geom_text(aes(label = Power), size = 3) +
  labs(x = "Soft threshold (power)", y = "Mean connectivity", title = "Mean connectivity") +
  theme_bw(base_size = 12)
ggsave(file.path(out, "figures", "D2_Fig1_softthreshold.pdf"), patchwork::wrap_plots(p1, p2, ncol = 2), width = 10, height = 4.2)
ggsave(file.path(out, "figures", "D2_Fig1_softthreshold.png"), patchwork::wrap_plots(p1, p2, ncol = 2), width = 10, height = 4.2, dpi = 300)

# 选择 power：fit R^2 >= 0.8 的最小 power；若无则取 R^2 最大值
signed_r2 <- -sign(sft_tab$slope) * sft_tab$SFT.R.sq
cand <- sft_tab$Power[signed_r2 >= 0.8]
softPower <- if (length(cand) > 0) min(cand) else sft_tab$Power[which.max(signed_r2)]
cat("soft threshold power =", softPower, " (max signed R2 =",
    round(max(signed_r2), 3), ")\n")

adj <- adjacency(datExpr, power = softPower, type = "signed")
TOM <- TOMsimilarity(adj, TOMType = "signed")
dissTOM <- 1 - TOM
geneTree <- hclust(as.dist(dissTOM), method = "average")

minModSize <- min(30, ncol(datExpr) - 5)
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModSize)
dynamicColors <- labels2colors(dynamicMods)
merged <- mergeCloseModules(datExpr, dynamicColors, cutHeight = 0.25, verbose = 3)
modColors <- merged$colors
modLabels <- labels2colors(modColors)
MEs <- merged$newMEs
cat("modules detected:", length(unique(modLabels)), "\n")
print(table(modLabels))

# 模块-性状相关
traits <- data.frame(
  StemnessScore = as.numeric(stem_score[rownames(datExpr)]),
  Lesion_binary = as.numeric(pheno[rownames(datExpr), "binary"] == "Lesion"),
  row.names = rownames(datExpr)
)
MEcor <- cor(MEs, traits, use = "p")
MEp   <- corPvalueStudent(MEcor, nrow(datExpr))
modtrait <- data.table(module = rownames(MEcor))
for (tname in colnames(MEcor)) {
  modtrait[[paste0("r_", tname)]]  <- round(MEcor[, tname], 4)
  modtrait[[paste0("p_", tname)]]  <- signif(MEp[, tname], 4)
}
modtrait[, moduleSize := as.integer(table(modLabels)[module])]
setcolorder(modtrait, c("module", "moduleSize"))
fwrite(modtrait, file.path(out, "results", "wgcna_module_trait_corr.tsv"), sep = "\t")
print(modtrait)

# 模块分配
mod_assign <- data.table(gene = colnames(datExpr), module = modLabels,
                         kME_Stemness = round(cor(datExpr, MEs[, 1], use = "p")[, 1], 4))
fwrite(mod_assign, file.path(out, "data", "wgcna_module_assignment.tsv"), sep = "\t")

# 模块-性状热图
textMatrix <- paste0(signif(MEcor, 2), "\n(", signif(MEp, 1), ")")
dim(textMatrix) <- dim(MEcor)
pdf(file.path(out, "figures", "D2_Fig2_module_trait_heatmap.pdf"), width = 8, height = 7)
par(mar = c(6, 9, 3, 3))
labeledHeatmap(Matrix = MEcor,
               xLabels = colnames(MEcor), yLabels = rownames(MEcor),
               ySymbols = rownames(MEcor), colorLabels = FALSE,
               colors = blueWhiteRed(50), textMatrix = textMatrix,
               setStdMargins = FALSE, cex.text = 0.7, zlim = c(-1, 1),
               main = "Module-trait relationships (GSE241715)")
dev.off()
png(file.path(out, "figures", "D2_Fig2_module_trait_heatmap.png"), width = 2400, height = 2100, res = 300)
par(mar = c(6, 9, 3, 3))
labeledHeatmap(Matrix = MEcor,
               xLabels = colnames(MEcor), yLabels = rownames(MEcor),
               ySymbols = rownames(MEcor), colorLabels = FALSE,
               colors = blueWhiteRed(50), textMatrix = textMatrix,
               setStdMargins = FALSE, cex.text = 0.7, zlim = c(-1, 1),
               main = "Module-trait relationships (GSE241715)")
dev.off()

# 模块树
pdf(file.path(out, "figures", "D2_Fig3_module_dendrogram.pdf"), width = 11, height = 6)
plotDendroAndColors(geneTree, modLabels, "Module",
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE,
                    guideHang = 0.05, main = "Gene dendrogram and module colours")
dev.off()
png(file.path(out, "figures", "D2_Fig3_module_dendrogram.png"), width = 3300, height = 1800, res = 300)
plotDendroAndColors(geneTree, modLabels, "Module",
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE,
                    guideHang = 0.05, main = "Gene dendrogram and module colours")
dev.off()

# -------------------------------------------------------------
# 5. StemnessScore 可视化
# -------------------------------------------------------------
stem_dt[, binary := factor(binary, levels = c("Normal","Lesion"))]
stem_dt[, group := factor(group, levels = c("Control", "Keratoconus", "Viral", "Bacterial_fungal"))]
pv <- ggplot(stem_dt, aes(group, StemnessScore, fill = group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.9) +
  labs(x = NULL, y = "StemnessScore (ssGSEA)", title = "Corneal limbal stemness score — GSE241715") +
  theme_bw(base_size = 12) + theme(legend.position = "none") +
  scale_fill_manual(values = c(Control = "#4C72B0", Keratoconus = "#DD8452",
                               Viral = "#55A868", Bacterial_fungal = "#C44E52"))
ggsave(file.path(out, "figures", "D2_Fig5_stemness_boxplot.pdf"), pv, width = 6, height = 4.2)
ggsave(file.path(out, "figures", "D2_Fig5_stemness_boxplot.png"), pv, width = 6, height = 4.2, dpi = 300)

# -------------------------------------------------------------
# 6. 验收判定
# -------------------------------------------------------------
sig <- modtrait[abs(r_StemnessScore) > 0.4 & p_StemnessScore < 0.05]
sig_les <- modtrait[abs(r_Lesion_binary) > 0.4 & p_Lesion_binary < 0.05]
verdict <- data.table(
  item = c("n_samples", "n_genes_pool", "soft_power", "n_modules",
           "modules_sig_StemnessScore", "modules_sig_Lesion", "pass_criterion"),
  value = c(nrow(datExpr), ncol(datExpr), softPower, length(unique(modLabels)),
            nrow(sig), nrow(sig_les),
            ">=1 module with |r|>0.4 & p<0.05 for both StemnessScore and Lesion")
)
fwrite(verdict, file.path(out, "results", "wgcna_summary.tsv"), sep = "\t")
print(verdict)

saveRDS(list(expr = datExpr, stem = stem_score, traits = traits, MEs = MEs,
             modLabels = modLabels, softPower = softPower, adj = adj, TOM = TOM,
             pheno = pheno, pool = pool_dt), file.path(out, "data", "D2_wgcna_objects.rds"))

cat("=== D2-02 done:", as.character(Sys.time()), "===\n")
sink(); close(log_con)