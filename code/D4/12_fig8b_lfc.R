# ============================================================
# D4-12 图 8B 重绘（GSE89537 面板 y 轴被 TMSB4X 拉伸，改为 log2FC 相对 day0）
# ============================================================
suppressPackageStartupMessages({ library(Biobase); library(data.table); library(ggplot2); library(patchwork) })
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
# 排除在 day0 表达 < 5 的基因（基线过低导致 lfc 不可解释）以及 TMSB4X（量级极端）
keep <- long3[, .(base = mean(base), mx = max(abs(lfc))), by = gene]
keep <- keep[base >= 5 & gene != "TMSB4X"]
long3b <- long3[gene %in% keep$gene]

pB <- ggplot(long3b, aes(timepoint, lfc, group = gene, colour = gene)) +
  geom_hline(yintercept = 0, linetype = 3, colour = "grey50") +
  geom_line(linewidth = 0.95) + geom_point(size = 2.6) +
  labs(title = "B  GSE89537 (mouse alkali-burn CNV): log2 fold change vs uninjured day 0",
       subtitle = "TMSB4X excluded (baseline ~11k, log2FC off-scale); baseline >= 5 retained; pooled group means (2 replicates/point)",
       x = NULL, y = "log2 FC vs Control_day0") +
  theme_bw(base_size = 10) +
  theme(legend.position = "right", legend.text = element_text(size = 7.5),
        legend.title = element_blank())

# 用已有 A 面板 PNG 拼不了，重跑 A 面板太贵 -> 直接单独输出 B 面板并与 A 用 patchwork 需 A 对象
# 方案：A 面板已在 09c 生成完整图；此处只替换 B 面板为独立文件，并在 09d 中合并
ggsave(file.path(out, "figures", "D4_Fig8B_mouse_lfc.pdf"), pB, width = 7.5, height = 4.6)
ggsave(file.path(out, "figures", "D4_Fig8B_mouse_lfc.png"), pB, width = 7.5, height = 4.6, dpi = 300)
print(long3b[, .(lfc_d5 = lfc[timepoint=="CNV_day5"], lfc_d15 = lfc[timepoint=="CNV_day15"]), by = gene])
cat("D4-12 done\n")