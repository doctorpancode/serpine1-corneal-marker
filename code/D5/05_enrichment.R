# ============================================================
# D5-05 功能富集：candidate pool + 核心基因高/低表达分层 GO/KEGG/Reactome
# 输出：results/D5_enrich_*.tsv
#       figures/D5_Fig10_enrichment.pdf/.png
# ============================================================
suppressPackageStartupMessages({
  library(Biobase); library(data.table); library(clusterProfiler)
  library(org.Hs.eg.db); library(msigdbr); library(ggplot2); library(patchwork)
})
out <- "D:/A生信分析/虚拟敲除/第5天执行结果"
DATA <- "D:/A生信分析/虚拟敲除/D1执行结果/D1_CNV_角膜新生血管_数据分析/data/processed"

pool <- fread("D:/A生信分析/虚拟敲除/第3天执行结果/data/candidate_pool.csv")$gene
cat("候选池:", length(pool), "\n")
ent <- suppressWarnings(bitr(pool, "SYMBOL", "ENTREZID", org.Hs.eg.db))
cat("映射成功:", nrow(ent), "\n")
univ <- suppressWarnings(bitr(rownames(readRDS(file.path(DATA, "GSE241715_expression.rds"))),
                              "SYMBOL", "ENTREZID", org.Hs.eg.db))$ENTREZID

# ---- GO BP / CC / MF + KEGG ----
go <- enrichGO(ent$ENTREZID, org.Hs.eg.db, ont = "ALL", universe = univ,
               pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
kegg <- tryCatch(enrichKEGG(ent$ENTREZID, organism = "hsa", universe = univ,
                            pvalueCutoff = 0.05, pAdjustMethod = "BH"), error = function(e) NULL)
cat("GO terms:", nrow(as.data.frame(go)), " KEGG:", if (is.null(kegg)) 0 else nrow(as.data.frame(kegg)), "\n")
fwrite(as.data.frame(go), file.path(out, "results", "D5_enrich_GO_all.tsv"), sep = "\t")
if (!is.null(kegg)) {
  kd <- as.data.frame(kegg)
  fwrite(kd, file.path(out, "results", "D5_enrich_KEGG.tsv"), sep = "\t")
}

# ---- Reactome（msigdbr C2:CP:REACTOME）----
msig <- msigdbr(species = "Homo sapiens", collection = "C2")
rt <- msig[grepl("^REACTOME_", msig$gs_name), ]
rtl <- split(rt$ncbi_gene, rt$gs_name)
rtl <- lapply(rtl, function(z) unique(as.character(z)))
rt_res <- tryCatch(enricher(ent$ENTREZID, TERM2GENE = data.frame(
  term = rep(names(rtl), lengths(rtl)), gene = unlist(rtl)), universe = univ,
  pvalueCutoff = 0.05, pAdjustMethod = "BH"), error = function(e) NULL)
if (!is.null(rt_res)) {
  rd <- as.data.frame(rt_res)
  setDT(rd)
  if ("ID" %in% names(rd)) setnames(rd, "ID", "term")
  rd[, term := gsub("^REACTOME_", "", term)]
  fwrite(rd, file.path(out, "results", "D5_enrich_Reactome.tsv"), sep = "\t")
  cat("Reactome terms:", nrow(rd), "\n")
}

# ---- 图10：GO 前 10 + KEGG 前 10 + Reactome 前 10 ----
gd <- as.data.frame(go)
gd <- gd[order(gd$p.adjust), ][1:min(30, nrow(gd)), ]
gd$Description <- ifelse(nchar(gd$Description) > 58,
                         paste0(substr(gd$Description, 1, 55), "..."), gd$Description)
gd$Description <- factor(gd$Description, levels = rev(unique(gd$Description)))
pG <- ggplot(gd, aes(-log10(p.adjust), Description, fill = ONTOLOGY)) +
  geom_col(width = 0.72) +
  scale_fill_brewer(palette = "Set2") +
  geom_vline(xintercept = -log10(0.05), linetype = 2, colour = "grey45") +
  labs(title = "GO enrichment of the 167-gene candidate pool", x = "-log10(adjusted P)",
       y = NULL, fill = "Ontology") +
  theme_bw(base_size = 9) + theme(legend.position = "top")
plts <- list(pG)
if (!is.null(kegg)) {
  kd <- as.data.frame(kegg); kd <- kd[order(kd$p.adjust), ][1:min(12, nrow(kd)), ]
  kd$Description <- factor(kd$Description, levels = rev(kd$Description))
  plts[[length(plts) + 1]] <- ggplot(kd, aes(-log10(p.adjust), Description)) +
    geom_col(width = 0.72, fill = "#8c96c6") +
    geom_vline(xintercept = -log10(0.05), linetype = 2, colour = "grey45") +
    labs(title = "KEGG pathway enrichment", x = "-log10(adjusted P)", y = NULL) +
    theme_bw(base_size = 9)
}
if (!is.null(rt_res)) {
  rd2 <- as.data.frame(rt_res)
  rd2 <- rd2[order(rd2$p.adjust), , drop = FALSE]
  rd2 <- head(rd2, 12)
  tcol <- if ("term" %in% names(rd2)) "term" else if ("ID" %in% names(rd2)) "ID" else names(rd2)[1]
  tt <- as.character(rd2[[tcol]])
  tt <- gsub("^REACTOME_", "", tt)
  tt <- ifelse(nchar(tt) > 52, paste0(substr(tt, 1, 49), "..."), tt)
  rd2$label <- factor(tt, levels = rev(tt))
  plts[[length(plts) + 1]] <- ggplot(rd2, aes(-log10(p.adjust), label)) +
    geom_col(width = 0.72, fill = "#66c2a5") +
    geom_vline(xintercept = -log10(0.05), linetype = 2, colour = "grey45") +
    labs(title = "Reactome pathway enrichment", x = "-log10(adjusted P)", y = NULL) +
    theme_bw(base_size = 9)
}
pE <- wrap_plots(plts, nrow = 1)
ggsave(file.path(out, "figures", "D5_Fig10_enrichment.pdf"), pE, width = 5.6 * length(plts), height = 5.6)
ggsave(file.path(out, "figures", "D5_Fig10_enrichment.png"), pE, width = 5.6 * length(plts), height = 5.6, dpi = 300)

cat("\n=== GO top10 ===\n"); print(gd[1:10, c("ONTOLOGY", "Description", "p.adjust", "Count")])
if (!is.null(kegg)) { cat("\n=== KEGG top10 ===\n"); print(as.data.frame(kegg)[1:10, c("Description", "p.adjust", "Count")]) }
if (!is.null(rt_res)) { cat("\n=== Reactome top10 ===\n"); print(head(as.data.frame(rt_res), 10)) }
cat("\nD5-05 done\n")