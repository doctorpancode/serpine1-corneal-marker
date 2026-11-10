# D3-02 候选池收口：与 limma 差异基因取交 + 支持度/模块分层
# 输入: data\candidate_pool.csv, D1 GSE241715_expression.rds
# 输出: data\candidate_pool.csv（覆盖为收口后版本，保留原始版为 candidate_pool_full.csv）
#       results\D3_limma_deg.tsv
#       results\D3_pool_tiering.tsv
#       figures\D3_Fig2_pool_volcano.pdf/.png
# =============================================================
suppressPackageStartupMessages({
  library(Biobase); library(limma); library(data.table); library(ggplot2); library(ggrepel)
})
d1  <- "D:/A生信分析/虚拟敲除/D1执行结果/D1_CNV_角膜新生血管_数据分析"
out <- "D:/A生信分析/虚拟敲除/第3天执行结果"

# 备份全量池
if (!file.exists(file.path(out, "data", "candidate_pool_full.csv")))
  file.copy(file.path(out, "data", "candidate_pool.csv"), file.path(out, "data", "candidate_pool_full.csv"))

eset <- readRDS(file.path(d1, "data/processed/GSE241715_expression.rds"))
expr <- exprs(eset); pd <- pData(eset)
grp  <- factor(ifelse(pd$group == "Control", "Normal", "Lesion"), levels = c("Normal", "Lesion"))

design <- model.matrix(~0 + grp); colnames(design) <- levels(grp)
fit <- lmFit(expr, design)
cm  <- makeContrasts(Lesion - Normal, levels = design)
fit2 <- eBayes(contrasts.fit(fit, cm))
tt <- topTable(fit2, number = Inf, sort.by = "P")
deg <- data.table(gene = rownames(tt), logFC = tt$logFC, AveExpr = tt$AveExpr,
                  t = tt$t, P = tt$P.Value, adjP = tt$adj.P.Val)
fwrite(deg, file.path(out, "results", "D3_limma_deg.tsv"), sep = "\t")
cat("limma: total", nrow(deg), " | adjP<0.05:", sum(deg$adjP < 0.05),
    " | |logFC|>0.5 & adjP<0.05:", sum(deg$adjP < 0.05 & abs(deg$logFC) > 0.5), "\n")

sig <- deg[adjP < 0.05 & abs(logFC) > 0.5]

# 读池
pool <- fread(file.path(out, "data", "candidate_pool_full.csv"))
pool <- merge(pool, deg[, .(gene, logFC, adjP, AveExpr)], by = "gene", all.x = TRUE)
cat("pool:", nrow(pool), " | 在 limma 中:", sum(!is.na(pool$logFC)), "\n")
pool[, is_DEG := as.integer(!is.na(adjP) & adjP < 0.05 & abs(logFC) > 0.5)]
cat("pool ∩ DEG:", sum(pool$is_DEG), "\n")

# 分层：T1 = 三主题全覆盖且是 DEG；T2 = 双主题且是 DEG；T3 = DEG；T4 = 单主题非DEG（保留备查）
pool[, n_theme := in_theme_A + in_theme_B + in_theme_C]
pool[, tier := fifelse(n_theme == 3 & is_DEG == 1, "T1_triple_DEG",
                fifelse(n_theme >= 2 & is_DEG == 1, "T2_multi_DEG",
                fifelse(is_DEG == 1, "T3_DEG",
                fifelse(n_theme == 3, "T4_triple_notDEG", "T5_other"))))]
print(pool[, .N, by = tier][order(tier)])
fwrite(pool[, .N, by = tier][order(tier)], file.path(out, "results", "D3_pool_tiering.tsv"), sep = "\t")

# 收口后的候选池：T1 ∪ T2 ∪ T3（即所有 DEG），另附 tier 与主题标记
final <- pool[is_DEG == 1]
final <- final[order(-n_theme, -n_source_sets)]
fwrite(final, file.path(out, "data", "candidate_pool.csv"), sep = ",")
cat("收口后候选池:", nrow(final), "\n")
cat("\n三个预备候选基因状态:\n")
print(pool[gene %in% c("SERPINE1", "THBS1", "LOX"),
           .(gene, n_theme, n_source_sets, tier, logFC = round(logFC, 3), adjP)])

# ---------- 火山图 ----------
deg[, cls := "ns"]
deg[adjP < 0.05 & abs(logFC) > 0.5, cls := "DEG"]
deg[, inpool := gene %in% final$gene]
lab <- deg[gene %in% c("SERPINE1", "THBS1", "LOX")]
p <- ggplot(deg, aes(logFC, -log10(P))) +
  geom_point(data = deg[cls == "ns"], colour = "grey80", size = 0.7, alpha = 0.6) +
  geom_point(data = deg[cls == "DEG" & inpool == FALSE], colour = "#4C72B0", size = 0.9, alpha = 0.7) +
  geom_point(data = deg[inpool == TRUE & cls == "DEG"], colour = "#C44E52", size = 1.3, alpha = 0.9) +
  geom_point(data = lab, colour = "black", size = 2.6, shape = 21, fill = "#F0E442") +
  geom_text_repel(data = lab, aes(label = gene), size = 3.4, fontface = "bold") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  labs(x = "log2 fold change (Lesion / Normal)", y = "-log10(P)",
       title = "GSE241715 differential expression",
       subtitle = "red = candidate pool genes; yellow labels = pre-specified candidates") +
  theme_bw(base_size = 12)
ggsave(file.path(out, "figures", "D3_Fig2_pool_volcano.pdf"), p, width = 7.4, height = 5.6)
ggsave(file.path(out, "figures", "D3_Fig2_pool_volcano.png"), p, width = 7.4, height = 5.6, dpi = 300)
cat("done\n")