# ============================================================
# D4-13 图8B 最终版：log2FC 面板 y 轴截断（GJA1/CAPG 极端值会压平其余基因）
# ============================================================
suppressPackageStartupMessages({ library(Biobase); library(data.table); library(ggplot2) })
out <- "D:/A生信分析/虚拟敲除/第4天执行结果"
genes <- c("TMSB4X","SERPINE1","S100A13","XIAP","CAPG","PRDX2","VDAC1","P4HB","SOD1","GJA1")
mmap <- c(TMSB4X="Tmsb4x", SERPINE1="Serpine1", S100A13="S100a13", XIAP="Xiap", CAPG="Capg",
          PRDX2="Prdx2", VDAC1="Vdac1", P4HB="P4hb", SOD1="Sod1", GJA1="Gja1")
e3 <- readRDS("D:/A生信分析/虚拟敲除/D1执行结果/D1_CNV_角膜新生血管_数据分析/data/processed/GSE89537_expression.rds")
ex3 <- exprs(e3); pd3 <- pData(e3)
mm <- unname(mmap[genes]); ok <- mm %in% rownames(ex3)
d3 <- as.data.table(t(ex3[mm[ok], , drop = FALSE]))
setnames(d3, mm[ok], genes[ok])
d3[, timepoint := factor(pd3$group, levels = c("Control_day0", "CNV_day5", "CNV_day15"))]
long3 <- melt(d3, id.vars = "timepoint", measure.vars = genes[ok],
              variable.name = "gene", value.name = "val")
base <- long3[timepoint == "Control_day0", .(base = mean(val)), by = gene]
long3 <- merge(long3, base, by = "gene")
long3[, lfc := val - base]
long3 <- long3[base >= 5 & gene != "TMSB4X"]

# 按 day5 的 |lfc| 排图例顺序
ord <- long3[timepoint == "CNV_day5"][order(-lfc)]$gene
long3[, gene := factor(gene, levels = rev(ord))]

pB <- ggplot(long3, aes(timepoint, lfc, group = gene, colour = gene)) +
  geom_hline(yintercept = 0, linetype = 3, colour = "grey50") +
  geom_line(linewidth = 1.0) + geom_point(size = 2.8) +
  scale_y_continuous(limits = c(-35, 60), oob = scales::squish,
                     breaks = seq(-30, 60, 15)) +
  labs(title = "B  GSE89537 (mouse alkali-burn CNV): log2 fold change vs uninjured day 0",
       subtitle = "pooled group means (2 biological replicates per point); TMSB4X excluded (baseline ~11k); CAPG day15 value clipped",
       x = NULL, y = "log2 FC vs Control_day0") +
  theme_bw(base_size = 10) +
  theme(legend.position = "right", legend.text = element_text(size = 7.5),
        legend.title = element_blank())
ggsave(file.path(out, "figures", "D4_Fig8B_mouse_lfc.pdf"), pB, width = 8.2, height = 4.8)
ggsave(file.path(out, "figures", "D4_Fig8B_mouse_lfc.png"), pB, width = 8.2, height = 4.8, dpi = 300)
cat("D4-13 done\n")