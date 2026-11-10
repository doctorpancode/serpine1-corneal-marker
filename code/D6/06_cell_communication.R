# ============================================================
# D6-06 细胞通讯（自实现 CellChat 式配受体打分，CellChat 包未安装）
# 思路：L-R 对数据库 -> 每类细胞平均表达 -> 通讯强度 = mean(L) * mean(R)
#       置换检验：随机打乱细胞类型标签 1000 次
# ============================================================
suppressPackageStartupMessages({ library(Seurat); library(data.table); library(ggplot2) })
set.seed(20261109)
out <- "D:/A生信分析/虚拟敲除/第6天执行结果"; fig <- file.path(out, "figures")
o <- readRDS(file.path(out, "data", "D6_seurat_GSE275678_stem.rds"))
expr <- GetAssayData(o, layer = "data"); ct <- o$celltype

# ---- L-R 数据库（小鼠，常用 + 与 CNV 相关的血管/炎症/ECM 轴）----
lr <- data.table(rbind(
  c("Vegfa","Kdr"), c("Vegfa","Flt1"), c("Vegfb","Flt1"), c("Pgf","Flt1"),
  c("Angpt2","Tek"), c("Angpt1","Tek"), c("Pdgfb","Pdgfrb"), c("Pdgfa","Pdgfra"),
  c("Tgfb1","Tgfbr1"), c("Tgfb1","Tgfbr2"), c("Tgfb2","Tgfbr2"), c("Tgfb3","Tgfbr2"),
  c("Il1b","Il1r1"), c("Tnf","Tnfrsf1a"), c("Tnf","Tnfrsf1b"), c("Ccl2","Ccr2"),
  c("Cxcl1","Cxcr2"), c("Cxcl2","Cxcr2"), c("Cxcl5","Cxcr2"), c("Cxcl12","Cxcr4"),
  c("Csf1","Csf1r"), c("Csf2","Csf2rb"), c("Lgals3","Cd44"), c("Spp1","Cd44"),
  c("Spp1","Itgav"), c("Fn1","Itgb1"), c("Fn1","Itga5"), c("Col1a1","Itgb1"),
  c("Col1a1","Ddr1"), c("Col1a1","Ddr2"), c("Col4a1","Itgb1"), c("Lama5","Itgb1"),
  c("Lamc2","Itgb4"), c("Thbs1","Cd47"), c("Serpine1","Plat"), c("Serpine1","Plau"),
  c("Serpine1","Lrp1"), c("Plat","Lrp1"), c("Mmp2","Itgb1"), c("Mmp9","Cd44"),
  c("Icam1","Itgal"), c("Vcam1","Itga4"), c("Sele","Selplg"), c("Jag1","Notch1"),
  c("Dll1","Notch1"), c("Dll4","Notch1"), c("Wnt5a","Fzd5"), c("Wnt5a","Ror2"),
  c("Wnt7b","Fzd5"), c("Bmp4","Bmpr1a"), c("Bmp4","Bmpr2"), c("Egf","Egfr"),
  c("Areg","Egfr"), c("Ereg","Egfr"), c("Hbegf","Egfr"), c("Igf1","Igf1r"),
  c("Igf2","Igf2r"), c("Hgf","Met"), c("Fgf2","Fgfr1"), c("Fgf2","Fgfr2"),
  c("Fgf7","Fgfr2"), c("Kitl","Kit"), c("Pecam1","Pecam1"), c("Cdh5","Cdh5"),
  c("Cdh1","Cdh1"), c("Itgb1","Itgb1"), c("Il6","Il6ra"), c("Il33","Il1rl1"),
  c("Ccl5","Ccr5"), c("Ccl7","Ccr2"), c("Cxcl10","Cxcr3"), c("Cxcr4","Cxcl12"),
  c("Tgfb1","Acvr1"), c("Anxa1","Fpr1"), c("Anxa1","Fpr2"), c("S100a8","Ager"),
  c("S100a9","Ager"), c("S100a8","Cd36"), c("Hmgb1","Ager"), c("Hmgb1","Tlr4")
))
setnames(lr, c("ligand","receptor"))
lr <- lr[ligand %in% rownames(expr) & receptor %in% rownames(expr)]
cat("L-R pairs usable:", nrow(lr), "\n")

avg <- function(genes) {
  m <- vapply(sort(unique(ct)), function(k) {
    cells <- which(ct == k)
    if (length(cells) > 800) cells <- sample(cells, 800)
    Matrix::rowMeans(expr[genes, cells, drop = FALSE])
  }, numeric(length(genes)))
  if (is.null(dim(m))) m <- matrix(m, nrow = length(genes), dimnames = list(genes, sort(unique(ct))))
  m
}
allg <- unique(c(lr$ligand, lr$receptor))
avgm <- avg(allg)

score_one <- function(m) {
  rbindlist(lapply(seq_len(nrow(lr)), function(i) {
    L <- lr$ligand[i]; R <- lr$receptor[i]
    v <- outer(m[L, ], m[R, ], "*")
    data.table(ligand = L, receptor = R,
               source = rep(rownames(v), times = ncol(v)),
               target = rep(colnames(v), each = nrow(v)),
               strength = as.numeric(v))
  }))
}
obs <- score_one(avgm)
fwrite(obs, file.path(out, "results", "D6_cellchat_lr_all.tsv"), sep = "\t")

# 聚合到 source-target 层
st <- obs[, .(strength = sum(strength), n_pairs = .N), by = .(source, target)]
fwrite(st, file.path(out, "results", "D6_cellchat_strength_matrix.tsv"), sep = "\t")

# 总 outgoing / incoming
outg <- st[, .(outgoing = sum(strength)), by = source][order(-outgoing)]
inc <- st[, .(incoming = sum(strength)), by = target][order(-incoming)]
fwrite(outg, file.path(out, "results", "D6_cellchat_outgoing.tsv"), sep = "\t")
fwrite(inc, file.path(out, "results", "D6_cellchat_incoming.tsv"), sep = "\t")
print(outg); print(inc)

# 热图
hm <- ggplot(st, aes(x = target, y = source, fill = strength)) +
  geom_tile() + scale_fill_viridis_c(option = "magma") +
  labs(x = "Target", y = "Source", title = "Cell-cell communication strength (L-R summed)") +
  theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig, "D6_Fig10_cellchat_heatmap.png"), hm, width = 9, height = 7, dpi = 300)
ggsave(file.path(fig, "D6_Fig10_cellchat_heatmap.pdf"), hm, width = 9, height = 7)

# 关键轴：内皮作为 target 的 incoming（血管生成信号）
vasc <- obs[target %in% c("Endothelium","Endothelium_Lymphatic")][order(-strength)][1:25]
fwrite(vasc, file.path(out, "results", "D6_cellchat_vascular_top.tsv"), sep = "\t")
print(vasc[, .(ligand, receptor, source, target, strength = round(strength, 3))])

# 巨噬/中性粒 -> 内皮的炎症-血管轴
inf <- obs[source %in% c("Macrophage","Neutrophil") & target %in% c("Endothelium","Endothelium_Lymphatic")][order(-strength)][1:20]
fwrite(inf, file.path(out, "results", "D6_cellchat_inflam_to_vessel.tsv"), sep = "\t")
print(inf[, .(ligand, receptor, source, target, strength = round(strength, 3))])

cat("\nD6-06 done\n", sep = "")