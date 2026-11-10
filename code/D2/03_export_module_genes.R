# =============================================================
# D2-03  WGCNA 显著性模块基因导出 + 与候选基因交叉核对
# 输入: data/D2_wgcna_objects.rds, data/wgcna_module_assignment.tsv,
#       results/wgcna_module_trait_corr.tsv
# 输出: data/wgcna_significant_module_genes.tsv
#       results/module_gene_counts.tsv
#       results/D2_candidate_gene_module_check.tsv
#       figures/D2_Fig6_sig_module_kME_strip.pdf/.png
# =============================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
out <- "D:/A生信分析/虚拟敲除/第2天执行结果"

objs <- readRDS(file.path(out, "data", "D2_wgcna_objects.rds"))
assign <- fread(file.path(out, "data", "wgcna_module_assignment.tsv"))
mt <- fread(file.path(out, "results", "wgcna_module_trait_corr.tsv"))
modsz <- table(factor(objs$modLabels))
mt[, moduleSize := as.integer(modsz[sub("^ME", "", module)])]
fwrite(mt, file.path(out, "results", "wgcna_module_trait_corr.tsv"), sep = "\t")

# 显著模块：与 StemnessScore 或 Lesion 至少一个 |r|>0.4 & p<0.05
sig <- mt[(abs(r_StemnessScore) > 0.4 & p_StemnessScore < 0.05) |
          (abs(r_Lesion_binary) > 0.4 & p_Lesion_binary < 0.05)]
sig[, module_short := sub("^ME", "", module)]
cat("significant modules:", nrow(sig), "\n"); print(sig)

genes_sig <- assign[module %in% sig$module_short]
genes_sig <- merge(genes_sig, sig[, .(module_short, module_full = module, r_StemnessScore,
                                     p_StemnessScore, r_Lesion_binary, p_Lesion_binary)],
                   by.x = "module", by.y = "module_short", all.x = TRUE)
genes_sig[, .kabs := abs(kME_Stemness)]
setorder(genes_sig, module, -.kabs)
genes_sig[, .kabs := NULL]
fwrite(genes_sig, file.path(out, "data", "wgcna_significant_module_genes.tsv"), sep = "\t")

cnt <- genes_sig[, .(n_genes = .N,
                     kME_abs_mean = round(mean(abs(kME_Stemness)), 4),
                     r_stemness = r_StemnessScore[1],
                     r_lesion = r_Lesion_binary[1]), by = module_full]
cnt[, .k := abs(r_stemness)]
setorder(cnt, -.k)
cnt[, .k := NULL]
fwrite(cnt, file.path(out, "results", "module_gene_counts.tsv"), sep = "\t")
print(cnt)

# 三个预备候选基因定位
cand <- c("SERPINE1", "THBS1", "LOX")
chk <- merge(data.table(gene = cand), assign, by = "gene", all.x = TRUE)
chk <- merge(chk, sig[, .(module_short, module_full = module, r_StemnessScore, r_Lesion_binary)],
             by.x = "module", by.y = "module_short", all.x = TRUE)
fwrite(chk, file.path(out, "results", "D2_candidate_gene_module_check.tsv"), sep = "\t")
print(chk)

# 图：显著模块的 kME 分布
p <- ggplot(genes_sig, aes(reorder(module_full, -abs(r_StemnessScore)), abs(kME_Stemness), fill = module_full)) +
  geom_violin(alpha = 0.55, scale = "width") +
  geom_boxplot(width = 0.12, outlier.size = 0.4, alpha = 0.9) +
  labs(x = NULL, y = "|kME| vs StemnessScore ME",
       title = "Connectivity of module genes to the stemness module eigengene") +
  theme_bw(base_size = 12) + theme(legend.position = "none")
ggsave(file.path(out, "figures", "D2_Fig6_sig_module_kME_strip.pdf"), p, width = 7, height = 4.2)
ggsave(file.path(out, "figures", "D2_Fig6_sig_module_kME_strip.png"), p, width = 7, height = 4.2, dpi = 300)
cat("done\n")