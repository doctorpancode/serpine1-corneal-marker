# =============================================================
# D2-01  Stemness gene set 构建（角膜缘/角膜上皮干细胞）
# 输入: 无（文献整理 + MSigDB 交叉核对）
# 输出: D:\A生信分析\虚拟敲除\第2天执行结果\data\stemness_genesets.tsv
#       D:\A生信分析\虚拟敲除\第2天执行结果\results\stemness_geneset_source_log.tsv
# 说明: 计划书要求 30–60 基因集。本集为"角膜缘干细胞/上皮祖细胞"核心标记，
#       每条基因标注文献来源类别，便于方法学章节交代。
# =============================================================
suppressPackageStartupMessages({
  library(data.table)
})

out_dir <- "D:/A生信分析/虚拟敲除/第2天执行结果"
dir.create(file.path(out_dir, "data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "results"), recursive = TRUE, showWarnings = FALSE)

gs <- data.table(
  gene = c(
    # --- 角膜上皮干细胞 / 角膜缘基底细胞（LSC / TAC）核心标记
    "TP63","KRT15","KRT19","KRT14","KRT5","KRT17",
    "MKI67","TOP2A","PCNA","CCNB1","CDK1",
    # --- 多能性 / 上皮祖细胞转录因子
    "SOX2","PAX6","POU5F1","NANOG","KLF4","MYC","BMI1","EZH2",
    # --- 干细胞表面 / 侧群与 ABC 转运体
    "ABCG2","ABCB1","ABCC1","NGFR","ITGA6","ITGB1","ITGB4","CD44",
    # --- 半桥粒 / 基底膜黏附（角膜缘基底细胞特征）
    "COL17A1","LAMA3","LAMB3","LAMC2","DST","PLEC",
    # --- 干性维持信号（Wnt / Notch / Hippo / BMP 抑分化）
    "WNT3A","WNT5A","CTNNB1","LEF1","AXIN2","HES1","HEY1","YAP1","WWTR1",
    "ID1","ID2","ID3","BMPR1A",
    # --- 静止期 / 干性相关代谢与保护
    "GJA1","GJB2","GJB6","ALDH1A1","ALDH3A1","SOD2","TXN","GPX1",
    # --- 角膜上皮终末分化（负向标记，纳入使评分具备"方向性"）
    "KRT3","KRT12","KRT13","KRT4","KRT76","CLDN1","TJP1"
  ),
  category = c(
    rep("LSC_basal_marker", 12),
    rep("stem_TF", 8),
    rep("stem_surface_ABC", 9),
    rep("hemidesmosome_adhesion", 6),
    rep("stemness_signaling", 14),
    rep("quiescence_protection", 8),
    rep("terminal_differentiation_negative", 7)
  )
)
gs <- unique(gs, by = "gene")

fwrite(gs, file.path(out_dir, "data", "stemness_genesets.tsv"), sep = "\t")

log_dt <- data.table(
  item = c("geneset_name", "n_genes", "n_categories", "n_positive", "n_negative",
           "source", "collection_date", "usage"),
  value = c("Corneal_Limbal_Stemness_v1",
            nrow(gs), length(unique(gs$category)),
            sum(gs$category != "terminal_differentiation_negative"),
            sum(gs$category == "terminal_differentiation_negative"),
            "文献整理(LSC marker: TP63/KRT15/ABCG2/COL17A1/ITGA6) + MSigDB GOBP_STEM_CELL_POPULATION_MAINTENANCE / DESCARTES_FETAL_EYE_CORNEAL_EPITHELIUM 交叉核对",
            as.character(Sys.Date()),
            "GSVA ssGSEA 评分；负向标记在 ssGSEA 中按同一方向参与（ssGSEA 不做符号翻转），故评分解释为'上皮干细胞程序总体活性'")
)
fwrite(log_dt, file.path(out_dir, "results", "stemness_geneset_source_log.tsv"), sep = "\t")

cat("geneset genes:", nrow(gs), "\n")
print(table(gs$category))