suppressPackageStartupMessages({ library(Biobase); library(data.table) })
es <- readRDS("D:/A生信分析/虚拟敲除/D1执行结果/D1_CNV_角膜新生血管_数据分析/data/processed/GSE241715_expression.rds")
x <- exprs(es); pd <- pData(es)
g <- c("TMSB4X","SERPINE1","S100A13","XIAP","CAPG","PRDX2","VDAC1","P4HB","SOD1","GJA1","THBS1","LOX")
d <- data.table(sample = colnames(x), group = pd$group, t(x[g, , drop = FALSE]))
setorder(d, group)
print(d)
cat("\n--- group means ---\n")
print(d[, lapply(.SD, mean), by = group, .SDcols = g])
cat("\n--- normal vs lesion wilcoxon p ---\n")
for (gg in g) {
  p <- wilcox.test(get(gg) ~ (group == "Control"), data = d)$p.value
  cat(sprintf("%-9s p=%.4g  meanCtrl=%.3f meanLesion=%.3f\n", gg, p,
              d[group == "Control", mean(get(gg))], d[group != "Control", mean(get(gg))]))
}