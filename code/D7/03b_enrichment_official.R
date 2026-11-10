# ============================================================
# D7-03b 官方虚拟敲除结果的富集（替换 RWR 版本）
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(clusterProfiler); library(org.Hs.eg.db)
  library(msigdbr); library(ggplot2); library(patchwork)
})
out <- "D:/A生信分析/虚拟敲除/第7天执行结果"
ko  <- fread(file.path(out, "results", "D7_knk_official_twoKO.tsv"))
cat("rows:", nrow(ko), " targets:", paste(unique(ko$target), collapse = ", "), "\n")

# 清单：每个敲除取全部显著基因（FDR<0.05，去掉目标自身）+ 前 50 名用来做背景
mklist <- function(g, n = 50, sig_only = FALSE) {
  d <- ko[target == g & gene != g][order(p.value)]
  if (sig_only) return(d[p.adj < 0.05]$gene)
  head(d$gene, n)
}
lists <- list(
  Serpine1_KO_sig = mklist("Serpine1", sig_only = TRUE),
  Tmsb4x_KO_sig   = mklist("Tmsb4x",   sig_only = TRUE),
  Serpine1_KO_top50 = mklist("Serpine1", 50),
  Tmsb4x_KO_top50   = mklist("Tmsb4x",   50),
  D3_pool_167 = fread("D:/A生信分析/虚拟敲除/第3天执行结果/data/candidate_pool.csv")$gene
)
for (nm in names(lists)) cat(nm, ":", length(lists[[nm]]), "genes\n")

m2h <- function(x) toupper(x)
all_res <- list()
for (nm in names(lists)) {
  g <- unique(m2h(lists[[nm]]))
  g <- g[!is.na(g) & nchar(g) > 1]
  if (length(g) < 5) { cat(nm, ": too few genes, skipped\n"); next }
  eg <- tryCatch(bitr(g, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
                 error = function(e) NULL)
  if (is.null(eg) || !nrow(eg)) next
  cat(nm, ": input", length(g), "mapped", nrow(eg), "\n")
  for (db in c("GO", "KEGG", "Reactome")) {
    r <- tryCatch({
      if (db == "GO") {
        enrichGO(eg$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE,
                 pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
      } else if (db == "KEGG") {
        enrichKEGG(eg$ENTREZID, organism = "hsa", pAdjustMethod = "BH",
                   pvalueCutoff = 0.05, qvalueCutoff = 0.2)
      } else {
        re <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
        tl <- split(re$gene_symbol, re$gs_name)
        enricher(g, TERM2GENE = data.frame(term = rep(names(tl), lengths(tl)),
                                           gene = unlist(tl, use.names = FALSE)),
                 pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
      }
    }, error = function(e) { cat("  ", db, "failed:", conditionMessage(e), "\n"); NULL })
    if (is.null(r) || !nrow(as.data.frame(r))) { cat("  ", db, ": 0 terms\n"); next }
    d <- as.data.table(as.data.frame(r)); d[, `:=`(list_name = nm, database = db)]
    all_res[[paste(nm, db)]] <- d
  }
}
enr <- rbindlist(all_res, fill = TRUE)
fwrite(enr, file.path(out, "results", "D7_enrichment_official_knk.tsv"), sep = "\t")
cat("\ntotal terms:", nrow(enr), "\n")
print(enr[, .N, by = .(list_name, database)])

# ---- 图：两个敲除的 GO/KEGG/Reactome top10 ----
plot_one <- function(d, ttl) {
  d <- d[order(p.adjust)][1:min(10, .N)]
  d[, Description := factor(substr(Description, 1, 58), levels = rev(substr(Description, 1, 58)))]
  ggplot(d, aes(x = Description, y = -log10(p.adjust), fill = -log10(p.adjust))) +
    geom_col(width = 0.72) + coord_flip() +
    scale_fill_gradient(low = "#9DC3E6", high = "#C00000") +
    labs(x = NULL, y = "-log10(FDR)", title = ttl) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none", plot.title = element_text(face = "bold", size = 10),
          axis.text.y = element_text(size = 8))
}
ps <- list()
for (nm in c("Serpine1_KO_sig", "Tmsb4x_KO_sig")) {
  for (db in c("GO", "KEGG", "Reactome")) {
    d <- enr[list_name == nm & database == db]
    if (nrow(d)) ps[[paste(nm, db)]] <- plot_one(d, sprintf("%s — %s", nm, db))
  }
}
if (length(ps)) {
  p <- wrap_plots(ps, ncol = 3)
  ggsave(file.path(out, "figures", "D7_Fig3b_enrichment_official.png"), p,
         width = 17, height = 4.0 * ceiling(length(ps) / 3), dpi = 300)
  ggsave(file.path(out, "figures", "D7_Fig3b_enrichment_official.pdf"), p,
         width = 17, height = 4.0 * ceiling(length(ps) / 3))
}
cat("\nD7-03b done\n")