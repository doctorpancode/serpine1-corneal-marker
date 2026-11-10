# =============================================================
# D2-05  StemnessScore 跨队列稳定性（GSE204791 人角膜 31 例）
# 计划书要求 D2 用 GSE204791 扩增正常对照并做跨队列稳定性检验
# 输入: D1 GSE204791_expression.rds + GSE204791_sample_sheet.tsv +
#       D1 results/GSE204791_probe_annotation.tsv
#       D2 data/stemness_genesets.tsv
# 输出: data/GSE204791_stemness_score.tsv
#       results/D2_crosscohort_stemness_test.tsv
#       figures/D2_Fig8_GSE204791_stemness_boxplot.pdf/.png
# =============================================================
suppressPackageStartupMessages({
  library(Biobase); library(GSVA); library(data.table); library(ggplot2)
})
d1  <- "D:/A生信分析/虚拟敲除/D1执行结果/D1_CNV_角膜新生血管_数据分析"
out <- "D:/A生信分析/虚拟敲除/第2天执行结果"

eset <- readRDS(file.path(d1, "data/processed/GSE204791_expression.rds"))
expr <- exprs(eset)
ss   <- fread(file.path(d1, "data/processed/GSE204791_sample_sheet.tsv"), sep = "\t")
ann  <- fread(file.path(d1, "results/GSE204791_probe_annotation.tsv"), sep = "\t")
cat("GSE204791:", nrow(expr), "probes x", ncol(expr), "samples\n")

# probe -> symbol（取行均值）
ann <- ann[!is.na(SYMBOL) & SYMBOL != ""]
m <- match(rownames(expr), ann$probe)
keep <- !is.na(m)
e2 <- expr[keep, , drop = FALSE]
sym <- ann$SYMBOL[m[keep]]
cat("mapped probes:", nrow(e2), "\n")
agg <- rowsum(e2, group = sym) / as.vector(table(sym)[sort(unique(sym))])
# rowsum 输出顺序为 sort(unique(sym))；用平均值
cnt <- as.vector(table(sym))
agg <- rowsum(e2, group = sym) / cnt[match(rownames(rowsum(e2, group = sym)), sort(unique(sym)))]
cat("genes after probe collapse:", nrow(agg), "\n")

gs <- fread(file.path(out, "data", "stemness_genesets.tsv"))
gl <- list(Corneal_Limbal_Stemness = intersect(unique(gs$gene), rownames(agg)))
cat("geneset present:", length(gl[[1]]), "\n")

par <- ssgseaParam(agg, gl, normalize = TRUE, minSize = 5, maxSize = 500)
sc <- suppressWarnings(gsva(par, verbose = FALSE))
dt <- data.table(sample = colnames(sc), StemnessScore = as.numeric(sc[1, ]))
dt <- merge(dt, ss[, .(gsm, compartment, group, age, gender, kc_stage)],
            by.x = "sample", by.y = "gsm", all.x = TRUE)
fwrite(dt, file.path(out, "data", "GSE204791_stemness_score.tsv"), sep = "\t")
print(dt[, .(n = .N, mean = round(mean(StemnessScore), 4), sd = round(sd(StemnessScore), 4)),
         by = .(compartment, group)][order(compartment, group)])

# 分compartment检验（KC vs Control）
res <- rbindlist(lapply(split(dt, dt$compartment), function(x) {
  w <- wilcox.test(x$StemnessScore[x$group == "KC"], x$StemnessScore[x$group == "Control"])
  data.table(compartment = x$compartment[1],
             n_control = sum(x$group == "Control"), n_kc = sum(x$group == "KC"),
             mean_control = round(mean(x$StemnessScore[x$group == "Control"]), 4),
             mean_kc = round(mean(x$StemnessScore[x$group == "KC"]), 4),
             diff = round(mean(x$StemnessScore[x$group == "KC"]) -
                            mean(x$StemnessScore[x$group == "Control"]), 4),
             p_value = signif(w$p.value, 4))
}))
print(res)
fwrite(res, file.path(out, "results", "D2_crosscohort_stemness_test.tsv"), sep = "\t")

# 与 GSE241715 方向一致性（同基因集 → 同方向即支持）
d2 <- fread(file.path(out, "data", "stemness_score.tsv"))
cat("\nGSE241715  mean Control =", round(mean(d2[group == "Control", StemnessScore]), 4),
    " mean Lesion =", round(mean(d2[group != "Control", StemnessScore]), 4), "\n")

dt[, grp := factor(paste(compartment, group, sep = "_"))]
p <- ggplot(dt, aes(grp, StemnessScore, fill = group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.12, size = 2) +
  labs(x = NULL, y = "StemnessScore (ssGSEA)",
       title = "Cross-cohort validation — GSE204791 (human cornea, 31 arrays)") +
  theme_bw(base_size = 12) +
  scale_fill_manual(values = c(Control = "#4C72B0", KC = "#C44E52")) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(out, "figures", "D2_Fig8_GSE204791_stemness_boxplot.pdf"), p, width = 6.2, height = 4.4)
ggsave(file.path(out, "figures", "D2_Fig8_GSE204791_stemness_boxplot.png"), p, width = 6.2, height = 4.4, dpi = 300)
cat("done\n")