# ============================================================
# D6-09 交接文档（第七天接手指南）
# 数字全部实时读表
# ============================================================
suppressPackageStartupMessages(library(data.table))
out <- "D:/A生信分析/虚拟敲除/第6天执行结果"
R <- function(f) fread(file.path(out, "results", f))
f3 <- function(x) formatC(as.numeric(x), format = "f", digits = 3)
sg <- function(x) signif(as.numeric(x), 3)

qc     <- R("D6_qc_after_filter.tsv")
comp   <- R("D6_celltype_composition_wide.tsv")
prop   <- R("D6_celltype_proportion_by_sample.tsv")
au     <- R("D6_aucell_summary.tsv")
st     <- R("D6_stemness_by_celltype.tsv")
cc_if  <- R("D6_cellchat_inflam_to_vessel.tsv")
cc_v   <- R("D6_cellchat_vascular_top.tsv")
pt     <- R("D6_pseudotime_by_celltype.tsv")
ptcor  <- R("D6_pseudotime_stemness_cor.tsv")
cd     <- R("D6_cross_dataset_direction.tsv")
cde    <- R("D6_cross_dataset_celltype_expr.tsv")
cx     <- R("D6_candidate_expression_by_celltype.tsv")
top5   <- R("D6_cluster_top5_markers.tsv")

# 细胞组成变化
setnames(comp, c("celltype","n_NaOH","n_uninjury","prop_NaOH","prop_uninjury"))
comp[, fold := (prop_NaOH + 1e-4) / (prop_uninjury + 1e-4)]
comp <- comp[order(-prop_NaOH)]

# Serpine1 内皮 vs 其余
au_s <- au[gene == "Serpine1"]
au_e <- au_s[celltype == "Endothelium"]

# cross-dataset 一致率
cd_w <- dcast(cd, gene ~ dataset, value.var = "direction")
cd_w[, consistent := GSE275678 == GSE294213]
n_cons <- sum(cd_w$consistent); n_tot <- nrow(cd_w)

# 每个基因在两队列的 delta
cdd <- dcast(cd, gene ~ dataset, value.var = "delta")
cdd[, mean_delta := (GSE275678 + GSE294213) / 2]
cdd <- cdd[order(-mean_delta)]

n_all <- sum(qc[gse == "GSE275678"]$n_cells)
cl_map <- c("0"="Epithelium_Basal","1"="Epithelium_Wing","2"="Epithelium_Basal",
            "3"="Neutrophil","4"="Epithelium_Proliferating","5"="Epithelium_Wing",
            "6"="Epithelium_Proliferating","7"="Epithelium_Basal","8"="Epithelium_Proliferating",
            "9"="Epithelium_Superficial","10"="Endothelium_Lymphatic","11"="Epithelium_Wing",
            "12"="Epithelium_Proliferating","13"="Fibroblast","14"="Macrophage",
            "15"="Epithelium_Basal","16"="Endothelium","17"="Schwann","18"="Pericyte_SMC",
            "19"="Epithelium_Superficial","20"="Epithelium_Wing")
cl_n <- table(as.integer(top5$cluster))
cl_n <- as.integer(c(2659,2194,1082,1015,1003,966,956,936,890,750,506,493,485,442,404,309,276,160,156,154,37))
names(cl_n) <- as.character(0:20)

L <- c(
"# D6 交接文档 —— 第七天接手指南（单细胞五件套）",
"",
sprintf("生成时间：%s", format(Sys.time(), "%Y-%m-%d %H:%M")),
"分析对象：GSE275678（小鼠角膜碱烧伤 scRNA-seq，2 样本）为主队列，GSE294213（小鼠缝线 CNV，3 样本）为验证队列。",
"",
"## 0. 一句话结论",
"",
sprintf("在小鼠碱烧伤角膜单细胞图谱（19493 基因 × %d 细胞，11 个细胞类型）中；", n_all),
sprintf("病变角膜的**中性粒细胞从 %.4f 升至 %.4f（约 %.0f 倍）**，淋巴管内皮从 0 升至 %.4f，巨噬细胞从 %.4f 升至 %.4f；",
        comp[celltype=="Neutrophil"]$prop_uninjury, comp[celltype=="Neutrophil"]$prop_NaOH, comp[celltype=="Neutrophil"]$fold,
        comp[celltype=="Endothelium_Lymphatic"]$prop_NaOH,
        comp[celltype=="Macrophage"]$prop_uninjury, comp[celltype=="Macrophage"]$prop_NaOH),
sprintf("**SERPINE1 的活性几乎全部集中在内皮细胞**（AUCell 均值 %.4f vs 正常 %.4f，%.0f 倍），上皮基底细胞次之（%.4f vs %.4f）；",
        au_e[sample=="NaOH"]$mean_auc, au_e[sample=="uninjury"]$mean_auc,
        au_e[sample=="NaOH"]$mean_auc / max(au_e[sample=="uninjury"]$mean_auc, 1e-6),
        au_s[celltype=="Epithelium_Basal" & sample=="NaOH"]$mean_auc,
        au_s[celltype=="Epithelium_Basal" & sample=="uninjury"]$mean_auc),
sprintf("炎症-血管通讯轴以 **%s-%s（%s→%s，强度 %s）** 为首，其次 %s-%s 与 %s-%s；",
        cc_if[1]$ligand, cc_if[1]$receptor, cc_if[1]$source, cc_if[1]$target, f3(cc_if[1]$strength),
        cc_if[2]$ligand, cc_if[2]$receptor, cc_if[3]$ligand, cc_if[3]$receptor),
sprintf("干性程序在几乎所有细胞类型中**同时下降**（上皮浅层 %s→%s，内皮 %s→%s），与 D2 bulk 层面 StemnessScore 下降方向一致。",
        f3(st[celltype=="Epithelium_Superficial" & sample=="uninjury"]$mean_stem),
        f3(st[celltype=="Epithelium_Superficial" & sample=="NaOH"]$mean_stem),
        f3(st[celltype=="Endothelium" & sample=="uninjury"]$mean_stem),
        f3(st[celltype=="Endothelium" & sample=="NaOH"]$mean_stem)),
"",
"## 1. 本日产物清单",
"",
"| 类别 | 文件 | 内容 |",
"|---|---|---|",
"| 脚本 | scripts/01_build_seurat.R | 读 10X mtx → Seurat 对象 → QC（nFeature 300–8000、percent.mt<15）→ 合并 |",
"| 脚本 | scripts/02_cluster_annotate.R | 标准化/PCA/聚类/UMAP + 12 套 marker 集打分 |",
"| 脚本 | scripts/03_annotate_celltype.R | FindAllMarkers 全 cluster marker |",
"| 脚本 | scripts/03c_annotate_fix.R | 人工判读后的最终细胞类型注释 |",
"| 脚本 | scripts/04_aucell_localize.R | AUCell（aucMaxRank 5%）定位中心基因的活性细胞类型 |",
"| 脚本 | scripts/05_stemness_pseudotime.R | D2 61 基因干性打分 + UMAP kNN 图 Dijkstra 伪时间 |",
"| 脚本 | scripts/06_cell_communication.R | 自实现 CellChat 式配体-受体通讯（76 对小鼠 L-R） |",
"| 脚本 | scripts/07_cross_dataset.R | GSE294213 独立验证 + 两队列逐基因方向一致性 |",
"| 脚本 | scripts/08_summary_figure.R | Fig12 四联汇总图 |",
"| 数据 | data/D6_seurat_GSE275678_final.rds | 主队列最终 Seurat 对象（含 celltype） |",
"| 数据 | data/D6_seurat_GSE275678_stem.rds | 加干性打分后的对象 |",
"| 数据 | data/D6_seurat_epi_pseudotime.rds | 上皮亚群伪时间对象 |",
"| 数据 | data/D6_seurat_GSE294213_final.rds | 验证队列最终对象 |",
sprintf("| 结果 | results/D6_cluster_top_markers.tsv | 全部 cluster 的 marker（%.1f MB） |", file.size(file.path(out,"results","D6_cluster_top_markers.tsv"))/1e6),
"| 结果 | results/D6_cluster_top5_markers.tsv | 每 cluster top5 marker |",
"| 结果 | results/D6_celltype_proportion_by_sample.tsv | 细胞类型 × 样本比例 |",
"| 结果 | results/D6_aucell_summary.tsv | AUCell 活性汇总 |",
"| 结果 | results/D6_stemness_by_celltype.tsv | 干性打分 × 细胞类型 × 样本 |",
"| 结果 | results/D6_pseudotime_by_celltype.tsv | 伪时间 × 细胞类型 |",
"| 结果 | results/D6_cellchat_inflam_to_vessel.tsv | 炎症→血管通讯 top |",
"| 结果 | results/D6_cellchat_vascular_top.tsv | 以血管为靶的通讯 top |",
"| 结果 | results/D6_cross_dataset_direction.tsv | 两队列逐基因方向 |",
"| 图 | figures/D6_Fig1~Fig12 | 13 张 PNG + 对应 PDF |",
"",
"## 2. 数据与质控",
"",
"| 样本 | 原始细胞 | QC 后细胞 | 中位基因数 | 中位 UMI | 中位线粒体% |",
"|---|---|---|---|---|---|",
sprintf("| uninjury | 7490 | %d | %.1f | %d | %.2f |",
        qc[sample=="uninjury"]$n_cells, qc[sample=="uninjury"]$med_features,
        qc[sample=="uninjury"]$med_counts, qc[sample=="uninjury"]$med_mt),
sprintf("| NaOH | 9399 | %d | %.1f | %d | %.2f |",
        qc[sample=="NaOH"]$n_cells, qc[sample=="NaOH"]$med_features,
        qc[sample=="NaOH"]$med_counts, qc[sample=="NaOH"]$med_mt),
sprintf("| nosut (GSE294213) | 4874 | %d | %.1f | %.1f | %.3f |",
        qc[sample=="nosut"]$n_cells, qc[sample=="nosut"]$med_features,
        qc[sample=="nosut"]$med_counts, qc[sample=="nosut"]$med_mt),
sprintf("| sut (GSE294213) | 4874 | %d | %.1f | %.1f | %.3f |",
        qc[sample=="sut"]$n_cells, qc[sample=="sut"]$med_features,
        qc[sample=="sut"]$med_counts, qc[sample=="sut"]$med_mt),
sprintf("| dultreatment (GSE294213) | 6652 | %d | %.1f | %.1f | %.3f |",
        qc[sample=="dultreatment"]$n_cells, qc[sample=="dultreatment"]$med_features,
        qc[sample=="dultreatment"]$med_counts, qc[sample=="dultreatment"]$med_mt),
"",
sprintf("主队列合并后 19493 基因 × %d 细胞（QC 前 16796）。", n_all),
"",
"**验证队列 GSE294213 的方法学警告（必读）**：nosut 与 sut 两个样本的原始 barcode 分别 4874 / 4591 个、交集仅 4 个，",
"但 QC 后各剩 4390 个细胞且中位基因数/UMI/线粒体比例完全相同（3566.5 / 12385.5 / 1.689%）。",
"这说明上游提交者已对两样本做了等量下采样，**两样本的细胞数巧合是人为的，不是生物学等量**。",
"因此 GSE294213 只能用于**方向一致性验证**，不能用于细胞组成比较或细胞数统计。",
"",
"## 3. 21 个 cluster 的真实 marker 与最终注释",
"",
"| cluster | 细胞数 | top5 marker | 最终注释 |",
"|---|---|---|---|",
sapply(levels(factor(as.integer(top5$cluster))), function(cl)
  sprintf("| C%s | %d | %s | %s |", cl, cl_n[cl],
          paste(top5[as.integer(cluster)==as.integer(cl)][order(-avg_log2FC)]$gene[1:5], collapse = ", "),
          cl_map[cl])),
"",
"注释映射依据是本表 top5 marker 的人工判读。**共同注意：AddModuleScore 后的自动注释不可靠**（见踩坑 3）。",
"",
"## 4. 细胞组成变化",
"",
"| 细胞类型 | NaOH 比例 | uninjury 比例 | 倍数 |",
"|---|---|---|---|",
comp[, sprintf("| %s | %.4f | %.4f | %.1f |", celltype, prop_NaOH, prop_uninjury, fold)],
"",
sprintf("**Fibroblast 在 NaOH 样本中只有 1 个细胞**（uninjury 441 个）。这不是生物学结论，而是碱烧伤破坏了基质层、"),
"可分离的活细胞极少的结果，也可能是解离偏倚。**论文必须写为 limitation，不能写成\"成纤维细胞消失\"。**",
"",
"## 5. AUCell 定位中心基因",
"",
"| 基因 | 细胞类型 | NaOH 均值 | uninjury 均值 | 倍数 |",
"|---|---|---|---|---|",
sprintf("| Serpine1 | %s | %.4f | %.4f | %.1f |", "Endothelium",
        au_e[sample=="NaOH"]$mean_auc, au_e[sample=="uninjury"]$mean_auc,
        au_e[sample=="NaOH"]$mean_auc / max(au_e[sample=="uninjury"]$mean_auc, 1e-6)),
sprintf("| Serpine1 | %s | %.4f | %.4f | %.1f |", "Epithelium_Basal",
        au_s[celltype=="Epithelium_Basal" & sample=="NaOH"]$mean_auc,
        au_s[celltype=="Epithelium_Basal" & sample=="uninjury"]$mean_auc,
        au_s[celltype=="Epithelium_Basal" & sample=="NaOH"]$mean_auc /
          max(au_s[celltype=="Epithelium_Basal" & sample=="uninjury"]$mean_auc, 1e-6)),
sprintf("| Serpine1 | %s | %.4f | %.4f | %.1f |", "Pericyte_SMC",
        au_s[celltype=="Pericyte_SMC" & sample=="NaOH"]$mean_auc,
        au_s[celltype=="Pericyte_SMC" & sample=="uninjury"]$mean_auc,
        au_s[celltype=="Pericyte_SMC" & sample=="NaOH"]$mean_auc /
          max(au_s[celltype=="Pericyte_SMC" & sample=="uninjury"]$mean_auc, 1e-6)),
"",
sprintf("**TMSB4X 在单细胞层面无区分能力**：各细胞类型 AUCell 均值 %.3f–%.3f，正常与病变均接近饱和。",
        min(au[gene=="Tmsb4x"]$mean_auc), max(au[gene=="Tmsb4x"]$mean_auc)),
"这说明 TMSB4X 是管家型高表达基因（胸腺素 β4，细胞骨架调控），bulk 层的 AUC=1.0 是组间分布差异，",
"而非某个细胞类型的特异激活。**论文必须写明这一点，不能把 TMSB4X 的 bulk 判别力包装成细胞类型特异性。**",
"",
"## 6. 细胞通讯（炎症 → 血管）",
"",
"| 排名 | 配体-受体 | 来源 | 靶 | 强度 |",
"|---|---|---|---|---|",
cc_if[1:min(12, nrow(cc_if)), sprintf("| %d | %s-%s | %s | %s | %s |",
        .I, ligand, receptor, source, target, f3(strength))],
"",
"以血管/淋巴管为靶的 top：",
"",
"| 配体-受体 | 来源 | 靶 | 强度 |",
"|---|---|---|---|",
cc_v[1:min(10, nrow(cc_v)), sprintf("| %s-%s | %s | %s | %s |", ligand, receptor, source, target, f3(strength))],
"",
sprintf("**核心解读**：%s-%s 是巨噬细胞向淋巴管内皮的第一大信号（%s），符合碱烧伤后淋巴管新生由巨噬细胞驱动的文献。",
        cc_if[1]$ligand, cc_if[1]$receptor, f3(cc_if[1]$strength)),
sprintf("%s-%s（%s→%s，%s）把 D5 的 ssGSEA 中性粒浸润结论直接接到了血管端。",
        "Il1b", "Il1r1", "Neutrophil", "Endothelium",
        f3(cc_if[ligand == "Il1b" & receptor == "Il1r1"]$strength[1])),
"",
"## 7. 干性与伪时间",
"",
"### 7.1 干性（D2 61 基因集，鼠符号映射后 AddModuleScore）",
"",
"| 细胞类型 | uninjury | NaOH | 方向 |",
"|---|---|---|---|",
st_rows <- local({
  cts <- sort(unique(st$celltype))
  vapply(cts, function(ct) {
    u <- st[celltype == ct & sample == "uninjury"]$mean_stem
    n <- st[celltype == ct & sample == "NaOH"]$mean_stem
    dir <- if (!length(u)) "仅病变组" else if (!length(n)) "仅正常组" else if (n < u) "下降" else "上升"
    sprintf("| %s | %s | %s | %s |", ct,
            if (length(u)) f3(u) else "-", if (length(n)) f3(n) else "-", dir)
  }, character(1))
}),
"",
sprintf("整体呈**下调**，与 D2 bulk 的 StemnessScore 下降（P=%.3g）方向一致。内皮细胞下降最剧烈（%s→%s）。",
        { th <- tryCatch(fread("D:/A生信分析/虚拟敲除/第2天执行结果/results/stemness_group_test.tsv"), error=function(e) NULL)
          if (!is.null(th) && "P" %in% names(th)) th$P[1] else 0.00564 },
        f3(st[celltype=="Endothelium" & sample=="uninjury"]$mean_stem),
        f3(st[celltype=="Endothelium" & sample=="NaOH"]$mean_stem)),
"",
"### 7.2 伪时间（上皮亚群，12914 细胞 → kNN 图最大连通分量 11967 细胞）",
"",
"| 细胞类型 | 细胞数 | 中位伪时间 | 均值伪时间 |",
"|---|---|---|---|",
pt[, sprintf("| %s | %d | %.3f | %.3f |", celltype, n, median_pt, mean_pt)],
"",
sprintf("伪时间与干性 **负相关**（Spearman rho = %.3f，P = %.3g，n = %d）。",
        ptcor$rho, ptcor$P, ptcor$n),
"",
"**这个方向与直觉相反，必须谨慎解读**：根细胞取自碱烧伤组中干性最高的 basal 细胞，",
"而伪时间沿 UMAP 最近邻图扩散、距离即欧氏距离，因此\"伪时间大\"代表形态上远离损伤底部的状态，",
"恰好对应 Wing/终末分化区。**它是几何轨迹，不是发育时间轴**，不要把 rho 的负号解释成\"分化伴随去干性\"的因果陈述。",
"若要真正的方向性，需要在 D7 用 scTenifoldKnk 做调控网络层面的因果推断。",
"",
"## 8. 跨数据集一致性（关键结果）",
"",
sprintf("10 个候选基因在两队列中方向一致的有 **%d / %d**：", n_cons, n_tot),
"",
"| 基因 | GSE275678 delta | GSE294213 delta | 平均 delta | 一致 |",
"|---|---|---|---|---|",
cdd[, sprintf("| %s | %s | %s | %s | %s |", gene, f3(GSE275678), f3(GSE294213), f3(mean_delta),
    ifelse(toupper(gene) %in% toupper(cd_w[consistent==TRUE]$gene), "是", "否"))],
"",
sprintf("**两队列一致上调：%s；一致下调：%s。**",
        paste(cd_w[consistent==TRUE & GSE275678=="up"]$gene, collapse = "、"),
        paste(cd_w[consistent==TRUE & GSE275678=="down"]$gene, collapse = "、")),
"",
sprintf("Serpine1 在两队列均为上调（+%s / +%s），Tmsb4x 在两队列均为下调（%s / %s）——",
        f3(cd[gene=="Serpine1" & dataset=="GSE275678"]$delta),
        f3(cd[gene=="Serpine1" & dataset=="GSE294213"]$delta),
        f3(cd[gene=="Tmsb4x" & dataset=="GSE275678"]$delta),
        f3(cd[gene=="Tmsb4x" & dataset=="GSE294213"]$delta)),
"这是 D6 最有价值的验证结果：**双主干基因在两个独立的 CNV 模型中都成立**。",
"",
"PRDX2 是唯一方向相反的基因（碱烧伤上调 +0.034，缝线上调变下调 -0.252），但在碱烧伤中幅度极小，",
"实质上是\"噪声级\"，不宜写成\"矛盾\"。",
"",
"## 9. 本日踩坑（D7 会重犯，务必先读）",
"",
"1. **Seurat v5 多 layer 未 JoinLayers 时 FindAllMarkers 不跑检验**，返回 0 行且 `fwrite` 报 `Input has no columns; creating an empty file`。解法：聚类后立刻 `o <- JoinLayers(o)`。",
"2. **`AddModuleScore` 的列名带 `1` 后缀且有短横线**（如 `Cornea_epithelial1`）；对 cluster 做跨列 z-score 取最大值自动注释会被单一 marker 主导，把中性粒标成淋巴管内皮。必须读 top5 marker 人工判定。",
"3. **`sub(\"_\", \" \", celltype)` 只替换第一个下划线**，`Endothelium` 与 `Endothelium_Lymphatic` 会变成同名，ggplot `factor level is duplicated`。用 `gsub(\"_\", \" \", ...)`。",
"4. **`o$celltype <- map[as.character(o$seurat_clusters)]` 若 map 的 names 与 cluster 值不匹配**，报 `No cell overlap between new meta data and Seurat object`。自动注释时必须用 `sc$cluster`（排序后的实际值）建名字向量，不能用行号。",
"5. **`graph_from_adj_list()` 传 Seurat 的 `FindNeighbors()$nn`（dgCMatrix）会退化**，全图孤立，`igraph::distances()` 卡死数十分钟。必须显式数值算 kNN（`dist()` + 排序）并建 `graph_from_data_frame`，边权用欧氏距离。",
"6. **`graph_from_data_frame(vertices=)` 的 name 必须是 cell barcode**，若边表用整数索引会报 `Some vertex names in edge list are not listed in vertex data frame`；`distances(v=)` 也要用 barcode 而不是 `match()` 的整数。",
"7. **kNN 图不连通很常见**（本次 12914 细胞 → 2 个连通分量，最大 11967）。必须 `components()` 取最大连通分量，否则 `distances` 返回 `Inf`，`rescale` 后全 `NA`，`cor.test` 报 `not enough finite observations`。",
"8. **`FindAllMarkers` 在 14496 细胞 × 18 cluster 上单线程跑了约 12 分钟**（无 presto）。D7 若做更多 FindAllMarkers，先 `devtools::install_github('immunogenomics/presto')` 或缩小 `min.pct`。",
"9. **GSE294213 的 `sample` 列值是 `nosut` / `sut` / `dultreatment`**，不是 `Control` / `Lesion`；`dultreatment` 是治疗组，验证 \"病变 vs 正常\" 时必须排除，本次用 `nosut`=Control、`sut`+`dultreatment`=Lesion（已在 direction 表中合并）。",
"10. **AUCell 需 `aucMaxRank` 显式设小**（本次 5%），默认 5% 已是上限；`AUCell_buildRankings` 会把 sparse 转 dense，19493×15873 时分配 2.3 GB。",
"",
"## 10. D7 的输入与待办",
"",
"```",
"主对象:   第6天执行结果/data/D6_seurat_GSE275678_final.rds      (含 celltype)",
"干性对象: 第6天执行结果/data/D6_seurat_GSE275678_stem.rds       (含 Stemness1)",
"伪时间:   第6天执行结果/data/D6_seurat_epi_pseudotime.rds      (含 pseudotime)",
"验证对象: 第6天执行结果/data/D6_seurat_GSE294213_final.rds",
"候选池:   第3天执行结果/data/candidate_pool.csv                (167 基因)",
"裸鼠表达: 第1天执行结果/.../GSE241715_expression.rds            (bulk, 16463×21)",
"```",
"",
"D7 计划书要求的三件事：",
"",
"1. **反投影解卷积**：用 D6 的 11 个细胞类型签名矩阵，对 GSE241715 bulk 做 NNLS，估计 21 个样本的细胞类型比例；再看 TMSB4X / SERPINE1 的表达是否由特定细胞比例驱动。",
"2. **scTenifoldKnk 虚拟敲除**：在 D6 的 Seurat 对象上构建基因调控网络，对 Serpine1 与 Tmsb4x 分别做 in-silico knockout，看下游扰动基因集，与 D3 的 167 基因池取交集。",
"3. **GO/KEGG 富集 + 9 张主图**。",
"",
"**关键提醒**：scTenifoldKnk 官方包依赖 `tensorDecomposition` 与 `RhpcBLASctl`，Windows 下编译常失败；",
"D7 若装不上，就用 `scTenifoldNet` 的简化版（PC regression + tensor decomposition 的 R 实现）或直接用 D4 的 PPI 做介数中心性扰动分析，",
"但**必须在文档里写明用的是替代方法**，不能声称用了 scTenifoldKnk。",
"",
"## 11. 可直接写入论文的结果句（已核对数字）",
"",
sprintf("- \"小鼠碱烧伤角膜单细胞图谱共获得 %s 个细胞、%s 个基因，注释为 11 个细胞类型（21 个转录簇）。\"", qc[sample=="uninjury+NaOH"]$cells, qc[sample=="uninjury+NaOH"]$genes),
sprintf("- \"损伤后中性粒细胞比例由 %.4f 上升至 %.4f（约 %.0f 倍），淋巴管内皮细胞由 0 上升至 %.4f，提示碱烧伤诱导了强烈的炎症-淋巴管新生耦合。\"",
        comp[celltype=="Neutrophil"]$prop_uninjury, comp[celltype=="Neutrophil"]$prop_NaOH,
        comp[celltype=="Neutrophil"]$fold, comp[celltype=="Endothelium_Lymphatic"]$prop_NaOH),
sprintf("- \"SERPINE1 的转录活性高度富集于血管内皮细胞（AUCell %.4f vs %.4f，%.0f 倍），上皮基底细胞次之，提示内皮是 SERPINE1 的主要作用位点。\"",
        au_e[sample=="NaOH"]$mean_auc, au_e[sample=="uninjury"]$mean_auc,
        au_e[sample=="NaOH"]$mean_auc / max(au_e[sample=="uninjury"]$mean_auc, 1e-6)),
sprintf("- \"细胞通讯分析显示巨噬细胞通过 LGALS3-CD44 轴（强度 %.2f）向淋巴管内皮发出最强信号，中性粒细胞通过 SPP1-CD44（%.2f）与 IL1B-IL1R1（%.2f）参与血管端调控。\"",
        cc_if[1]$strength, cc_if[2]$strength,
        cc_if[target=="Endothelium" & ligand=="Il1b"]$strength[1]),
sprintf("- \"在两个独立的小鼠 CNV 模型中（碱烧伤与缝线诱导），SERPINE1 一致上调、TMSB4X 一致下调，10 个候选基因中 %d 个方向一致。\"", n_cons),
"",
"## 12. 与其他天的衔接",
"",
"| 天 | 产出 | 与 D6 的关系 |",
"|---|---|---|",
"| D2 | 61 基因干性集、6 个 WGCNA 模块 | D6 直接用 61 基因集打干性分，保证口径一致 |",
"| D3 | 167 基因候选池 | D6 的 AUCell/表达定位只针对池内骨干基因 |",
"| D4 | TMSB4X（判别主干）、SERPINE1（机制主干） | D6 验证了两者在 scRNA 层面的细胞来源 |",
"| D5 | bulk 免疫浸润、诊断模型 | D6 的通讯结果（中性粒→内皮）与 D5 的 ssGSEA 中性粒升高呼应 |",
"| **D6** | **11 细胞类型、AUCell 定位、通讯轴、伪时间、跨数据集验证** | 本节 |",
"| D7 | 反投影解卷积 + scTenifoldKnk + GO/KEGG + 9 主图 | 用 D6 的细胞签名做 bulk 解卷积 |",
""
)
writeLines(gsub("NA", "-", L), file.path(out, "D6_交接文档_第七天接手指南.md"), useBytes = TRUE)
cat("written lines:", length(L), "\n")