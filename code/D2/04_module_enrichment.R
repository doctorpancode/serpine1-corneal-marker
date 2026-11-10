# =============================================================
# D2-04  显著模块 GO/KEGG 富集（为 D3 候选池提供"模块基因"视角）
# 输入: data/wgcna_significant_module_genes.tsv
# 输出: results/D2_module_enrichment_GO.tsv, results/D2_module_enrichment_KEGG.tsv
#       figures/D2_Fig7_module_enrichment.pdf/.png
# =============================================================
suppressPackageStartupMessages({
  library(data.table); library(clusterProfiler); library(org.Hs.eg.db)
  library(ggplot2); library(enrichplot)
})
out <- "D:/A生信分析/虚拟敲除/第2天执行结果"
set.seed(20261109)

gs <- fread(file.path(out, "data", "wgcna_significant_module_genes.tsv"))
cat("modules:", paste(unique(gs$module_full), collapse = ", "), "\n")

universe <- keys(org.Hs.eg.db, keytype = "SYMBOL")
res_list_go <- list(); res_list_kegg <- list()

for (m in unique(gs$module_full)) {
  g <- unique(gs[module_full == m, gene])
  eg <- bitr(g, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  cat("\n##", m, "n =", length(g), "mapped =", nrow(eg), "\n")
  if (nrow(eg) < 10) next
  go <- enrichGO(gene = eg$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                 ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                 readable = TRUE, universe = NULL)
  if (!is.null(go) && nrow(as.data.frame(go)) > 0) {
    d <- as.data.frame(go); d$module <- m
    res_list_go[[m]] <- as.data.table(d)
  }
  kg <- tryCatch(enrichKEGG(gene = eg$ENTREZID, organism = "hsa",
                            pvalueCutoff = 0.05, qvalueCutoff = 0.2),
                 error = function(e) NULL)
  if (!is.null(kg) && nrow(as.data.frame(kg)) > 0) {
    d <- as.data.frame(kg)
    if (!is.null(d$geneID)) d$geneID_symbol <- sapply(strsplit(d$geneID, "/"), function(x) {
      m2 <- eg$SYMBOL[match(x, eg$ENTREZID)]; paste(na.omit(m2), collapse = "/") })
    d$module <- m
    res_list_kegg[[m]] <- as.data.table(d)
  }
}
go_all <- rbindlist(res_list_go, fill = TRUE)
kg_all <- rbindlist(res_list_kegg, fill = TRUE)
fwrite(go_all, file.path(out, "results", "D2_module_enrichment_GO.tsv"), sep = "\t")
fwrite(kg_all, file.path(out, "results", "D2_module_enrichment_KEGG.tsv"), sep = "\t")
cat("\nGO rows:", nrow(go_all), " KEGG rows:", nrow(kg_all), "\n")

# 图：各模块 top5 GO BP
if (nrow(go_all) > 0) {
  top <- go_all[, head(.SD[order(pvalue)], 5), by = module]
  top[, term_short := substr(Description, 1, 55)]
  top[, score := -log10(pvalue)]
  p <- ggplot(top, aes(reorder(term_short, score), score, fill = module)) +
    geom_col() + coord_flip() + facet_wrap(~module, scales = "free_y", ncol = 2) +
    labs(x = NULL, y = "-log10(P)", title = "Top GO BP terms per WGCNA module") +
    theme_bw(base_size = 10) + theme(legend.position = "none")
  ggsave(file.path(out, "figures", "D2_Fig7_module_enrichment.pdf"), p, width = 11, height = 8)
  ggsave(file.path(out, "figures", "D2_Fig7_module_enrichment.png"), p, width = 11, height = 8, dpi = 300)
}
cat("done\n")