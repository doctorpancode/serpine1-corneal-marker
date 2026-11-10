# ============================================================
# D4-10 人-鼠同源基因映射（GSE89537 行名是小鼠符号/Ensembl 风格，如 Serpine1）
# 用 NCBI gene_info 的 Mammalia 子集？不可用 -> 转用 msigdbr / homologene 不行
# 方案: 用 R 包 org.Hs.eg.db + org.Mm.eg.db? 未装 org.Mm.eg.db
#       改用 babelgene / homologene? 需联网安装
# 最终方案: 直接按大小写转换匹配（小鼠符号 = 首字母大写 + 其余小写），
#           并对 GSE89537 的 "A,B" 双符号行拆分后逐一匹配
# ============================================================
suppressPackageStartupMessages({ library(Biobase); library(data.table) })
e3 <- readRDS("D:/A生信分析/虚拟敲除/D1执行结果/D1_CNV_角膜新生血管_数据分析/data/processed/GSE89537_expression.rds")
rn <- rownames(e3)
cat("矩阵行数:", length(rn), "\n")
cat("首字母大写比例:", round(mean(grepl("^[A-Z][a-z0-9]+$", rn)), 4), "\n")
cat("含逗号行数:", sum(grepl(",", rn)), "\n")
g <- c("TMSB4X","SERPINE1","S100A13","XIAP","CAPG","PRDX2","VDAC1","P4HB","SOD1","GJA1")
mm <- c("Tmsb4x","Serpine1","S100a13","Xiap","Capg","Prdx2","Vdac1","P4hb","Sod1","Gja1")
hit <- mm %in% rn
cat("小鼠符号命中:", sum(hit), "/", length(mm), "\n")
print(mm[hit])
print(setdiff(mm, rn))
# 也检查是否是 "0610005C13Rik,Serpine1" 这类合并行
for (x in mm) {
  cand <- rn[grepl(paste0("(^|,)", x, "($|,)"), rn)]
  if (length(cand)) cat(x, "->", paste(cand, collapse = " | "), "\n")
}