# ============================================================
# D7-03 GO/KEGG/Reactome 富集：三条基因清单分别富集
#   (1) Serpine1 KO top100（RWR 影响）
#   (2) Tmsb4x  KO top100
#   (3) D3 167 候选池（对照，确认与 D3 一致）
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(clusterProfiler); library(org.Hs.eg.db); library(ggplot2); library(msigdbr)
})
out <- "D:/A生信分析/虚拟敲除/第7天执行结果"
fig <- file.path(out, "figures")
ko  <- fread(file.path(out, "results", "D7_knockout_rwr.tsv"))

m2h <- function(x) toupper(x)   # 小鼠符号 -> 人类大写

lists <- list(
  Serpine1_KO_top100 = head(ko[target == "Serpine1" & gene != "Serpine1"][order(-abs(delta))]$gene, 100),
  Tmsb4x_KO_top100   = head(ko[target == "Tmsb4x"  & gene != "Tmsb4x" ][order(-abs(delta))]$gene, 100)
)
pool <- fread("D:/A生信分析/虚拟敲除/第3天执行结果/data/candidate_pool.csv")
lists$D3_pool_167 <- pool$gene

all_res <- list()
for (nm in names(lists)) {
  g <- unique(m2h(lists[[nm]]))
  g <- g[!is.na(g) & nchar(g) > 1]
  eg <- tryCatch(bitr(g, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
                 error = function(e) NULL)
  if (is.null(eg) || !nrow(eg)) { cat(nm, ": no mapped genes\n"); next }
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
        t2g <- data.frame(term = rep(names(tl), lengths(tl)), gene = unlist(tl, use.names = FALSE))
        enricher(g, TERM2GENE = t2g, pAdjustMethod = "BH",
                 pvalueCutoff = 0.05, qvalueCutoff = 0.2)
      }
    }, error = function(e) { cat(nm, db, "failed:", conditionMessage(e), "\n"); NULL })
    if (is.null(r) || !nrow(as.data.frame(r))) next
    d <- as.data.table(as.data.frame(r))
    d[, `:=`(list_name = nm, database = db)]
    all_res[[paste(nm, db)]] <- d
  }
}
enr <- rbindlist(all_res, fill = TRUE)
fwrite(enr, file.path(out, "results", "D7_enrichment_all.tsv"), sep = "\t")
cat("\ntotal enriched terms:", nrow(enr), "\n")
print(enr[, .N, by = .(list_name, database)])

# ---- 图：每张清单 top10 ----
plot_one <- function(d, ttl) {
  d <- d[order(p.adjust)][1:min(10, .N)]
  d[, Description := factor(Description, levels = rev(Description))]
  ggplot(d, aes(x = Description, y = -log10(p.adjust), fill = -log10(p.adjust))) +
    geom_col(width = 0.72) + coord_flip() +
    scale_fill_gradient(low = "#9DC3E6", high = "#C00000") +
    labs(x = NULL, y = "-log10(FDR)", title = ttl) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", hjust = 0),
          axis.text.y = element_text(size = 8))
}
ps <- list()
for (nm in unique(enr$list_name)) {
  for (db in unique(enr[list_name == nm]$database)) {
    d <- enr[list_name == nm & database == db]
    if (nrow(d)) ps[[paste(nm, db)]] <- plot_one(d, sprintf("%s — %s", nm, db))
  }
}
suppressPackageStartupMessages(library(patchwork))
if (length(ps) >= 2) {
  p <- wrap_plots(ps, ncol = min(3, length(ps)))
  h <- 4.2 * ceiling(length(ps) / min(3, length(ps)))
  ggsave(file.path(fig, "D7_Fig3_enrichment.png"), p, width = 16, height = max(6, h), dpi = 300)
  ggsave(file.path(fig, "D7_Fig3_enrichment.pdf"), p, width = 16, height = max(6, h))
}
cat("\nD7-03 done\n")