# ============================================================
# D7-03c 官方敲除富集图（只用官方 scTenifoldKnk 输出的两个敲除）
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})
out <- "D:/A生信分析/虚拟敲除/第7天执行结果"
e <- fread(file.path(out, "results", "D7_enrichment_official_knk.tsv"))
cat("terms:", nrow(e), "\n")
print(e[, .N, by = .(list_name, database)])

e <- e[list_name %in% c("Serpine1_KO_top50", "Tmsb4x_KO_top50")]
e[, nratio := as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio))]
e[, db2 := factor(database, levels = c("GO", "KEGG", "Reactome"))]

plot_one <- function(d, ttl) {
  d <- d[order(p.adjust)][1:min(10, .N)]
  lab <- ifelse(nchar(d$Description) > 60, paste0(substr(d$Description, 1, 57), "..."), d$Description)
  d[, lab := factor(lab, levels = rev(unique(lab)))]
  ggplot(d, aes(x = lab, y = -log10(p.adjust), fill = nratio)) +
    geom_col(width = 0.72) + coord_flip() +
    scale_fill_gradient(low = "#BDD7EE", high = "#C00000", name = "Gene\nratio") +
    labs(x = NULL, y = "-log10(FDR)", title = ttl) +
    theme_bw(base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 9.5),
          axis.text.y = element_text(size = 7.2),
          legend.key.size = unit(0.28, "cm"),
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 6.5),
          plot.margin = margin(4, 6, 4, 4))
}

ps <- list()
for (ko in c("Serpine1_KO_top50", "Tmsb4x_KO_top50")) {
  for (db in c("GO", "KEGG", "Reactome")) {
    d <- e[list_name == ko & database == db]
    tag <- ifelse(ko == "Serpine1_KO_top50", "A", "B")
    if (nrow(d) == 0) {
      ps[[paste(ko, db)]] <- ggplot() + theme_void() +
        annotate("text", x = 0, y = 0, label = "no term", size = 3, colour = "grey50")
      next
    }
    ps[[paste(ko, db)]] <- plot_one(d, sprintf("%s  %s — %s", tag, sub("_KO_top50", "", ko), db))
  }
}
p <- wrap_plots(ps, ncol = 3)
ggsave(file.path(out, "figures", "D7_Fig3b_enrichment_official.png"), p,
       width = 16.5, height = 7.4, dpi = 300, limitsize = FALSE)
ggsave(file.path(out, "figures", "D7_Fig3b_enrichment_official.pdf"), p,
       width = 16.5, height = 7.4, limitsize = FALSE)
cat("written\n")

# 打印文本版 top 供写文档
for (ko in c("Serpine1_KO_top50", "Tmsb4x_KO_top50")) {
  cat("\n########", ko, "########\n")
  for (db in c("GO", "KEGG", "Reactome")) {
    d <- e[list_name == ko & database == db][order(p.adjust)][1:min(6, .N)]
    if (!nrow(d)) next
    cat("\n--", db, "--\n")
    for (i in seq_len(nrow(d))) cat(sprintf("  %-62s %6s  FDR=%.3g\n", substr(d$Description[i], 1, 62), d$GeneRatio[i], d$p.adjust[i]))
  }
}
cat("\nD7-03c done\n")