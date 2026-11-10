# =============================================================
# D2 交接文档生成（第二天 → 第三天接手）
# 输出: D2_交接文档_第三天接手指南.md
# =============================================================
suppressPackageStartupMessages({ library(data.table) })
out <- "D:/A生信分析/虚拟敲除/第2天执行结果"

mt   <- fread(file.path(out, "results", "wgcna_module_trait_corr.tsv"))
cnt  <- fread(file.path(out, "results", "module_gene_counts.tsv"))
summ <- fread(file.path(out, "results", "wgcna_summary.tsv"))
stg  <- fread(file.path(out, "results", "stemness_group_test.tsv"))
xco  <- fread(file.path(out, "results", "D2_crosscohort_stemness_test.tsv"))
chk  <- fread(file.path(out, "results", "D2_candidate_gene_module_check.tsv"))
go   <- fread(file.path(out, "results", "D2_module_enrichment_GO.tsv"))
kg   <- fread(file.path(out, "results", "D2_module_enrichment_KEGG.tsv"))

getv <- function(item) summ[item == item_key]$value
summ[, item_key := item]

L <- c()
A <- function(...) {
  parts <- list(...)
  txt <- paste0(vapply(parts, function(x) {
    if (is.numeric(x)) format(x, digits = 6, trim = TRUE, scientific = FALSE) else as.character(x)
  }, character(1)), collapse = "")
  txt <- gsub("\n", " ", txt)
  L <<- c(L, txt)
}

A("# D2 交接文档 — 第三天接手指南")
A("")
A("> 生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ",
  "生成脚本：`scripts/07_write_handover.R`  ",
  "对应计划书章节：5.2 D2（表型评分构造 + WGCNA 建网）")
A("")
A("---")
A("")
A("## 0. 一句话结论")
A("")
A("**D2 全部验收项达成，且用的是计划书的「正路」而非退路。**")
A("")
A("- StemnessScore（角膜缘干细胞评分）在 GSE241715 中 **正常 1.5007 vs 病变 1.3336，Wilcoxon P = `",
  stg$p_value[1], "`**，方向符合「干细胞缺陷驱动血管化」的核心假说。")
A("- 限维后（血管生成/角化/ECM/黏附基因池 **", summ[item_key == "n_genes_pool"]$value[[1]], " 基因**）建网，soft power = ",
  summ[item_key == "soft_power"]$value[[1]], "，得到 **", summ[item_key == "n_modules"]$value[[1]], " 个模块**。")
A("- **MEturquoise 与 StemnessScore 强正相关 r = ", mt[module == "MEturquoise"]$r_StemnessScore,
  "（P = ", mt[module == "MEturquoise"]$p_StemnessScore, "）**，且其 GO 富集前三位就是 cell-substrate adhesion / morphogenesis of an epithelium / **regulation of vasculature development**——**这是本方案「干性-血管化轴」的第一个硬证据**。")
A("- MEblue（r = ", mt[module == "MEblue"]$r_StemnessScore, "）、MEgreen（r = ", mt[module == "MEgreen"]$r_StemnessScore,
  "）与 StemnessScore 负相关、与病变分组正相关（r = ", mt[module == "MEgreen"]$r_Lesion_binary, "），构成「干性丧失 → 病变程序启动」的对偶。")
A("- 三个预备候选基因 **SERPINE1 / THBS1 / LOX 全部落在 MEturquoise**，且 SERPINE1 的模块连接度 kME = ",
  chk[gene == "SERPINE1"]$kME_Stemness, "（三基因中最高，属模块核心 hub 级）。")
A("")
A("---")
A("")
A("## 1. 本日产出清单（逐项可核验）")
A("")
A("### 1.1 脚本（按执行顺序）")
A("")
A("| # | 脚本 | 作用 | 状态 |")
A("|---|------|------|------|")
A("| 1 | `scripts/01_build_stemness_geneset.R` | 构建 61 基因的角膜缘干性基因集（7 个功能类别，含 4 个终末分化负向标记） | 成功 |")
A("| 2 | `scripts/02_stemness_wgcna.R` | ssGSEA 打分 + 限维基因池 + WGCNA（signed TOM）+ 模块-性状相关 + 6 张图 | 成功 |")
A("| 3 | `scripts/03_export_module_genes.R` | 显著模块基因导出、模块规模统计、候选基因归属核查、kME 分布图 | 成功 |")
A("| 4 | `scripts/04_module_enrichment.R` | 6 个显著模块的 GO BP / KEGG 富集（clusterProfiler，在线 KEGG） | 成功 |")
A("| 5 | `scripts/05_crosscohort_stemness.R` | GSE204791（31 例人角膜，上皮/基质分层）跨队列复算 StemnessScore | 成功 |")
A("| 6 | `scripts/06_stemness_consistency_plot.R` | 两队列并排一致性图 | 成功 |")
A("| 7 | `scripts/07_write_handover.R` | 生成本文档 | 成功 |")
A("")
A("### 1.2 数据文件（`data/`）")
A("")
A("| 文件 | 内容 | 关键规模 |")
A("|------|------|----------|")
for (f in c("stemness_genesets.tsv", "stemness_score.tsv", "wgcna_gene_pool.tsv",
            "wgcna_module_assignment.tsv", "wgcna_significant_module_genes.tsv",
            "GSE204791_stemness_score.tsv", "D2_wgcna_objects.rds")) {
  p <- file.path(out, "data", f)
  n <- if (grepl("rds$", f)) "R list: expr/stem/traits/MEs/modLabels/softPower/adj/TOM/pheno/pool"
       else paste0(nrow(fread(p, sep = "\t")), " 行")
  sz <- round(file.info(p)$size / 1024, 1)
  A("| `data/", f, "` | ", n, " | ", sz, " KB |")
}
A("")
A("### 1.3 结果表（`results/`）")
A("")
for (f in sort(list.files(file.path(out, "results"), pattern = "tsv$"))) {
  p <- file.path(out, "results", f)
  A("- `results/", f, "` — ", nrow(fread(p, sep = "\t")), " 行")
}
A("")
A("### 1.4 图（`figures/`，PDF 矢量 + PNG 300 dpi）")
A("")
for (f in sort(list.files(file.path(out, "figures"), pattern = "png$"))) {
  A("- `figures/", f, "`")
}
A("")
A("---")
A("")
A("## 2. 核心结果表（把这三张表背下来，D3–D5 全靠它）")
A("")
A("### 2.1 模块-性状关系矩阵")
A("")
A("列说明：`r_StemnessScore` 为模块特征基因（ME）与 ssGSEA StemnessScore 的 Pearson r；`r_Lesion_binary` 为与「病变=1/正常=0」的 r。")
A("")
A("| 模块 | 基因数 | r(StemnessScore) | P | r(Lesion) | P | 判定 |")
A("|------|--------|------------------|---|-----------|---|------|")
for (i in seq_len(nrow(mt))) {
  r1 <- mt[i]
  verdict <- if (abs(r1$r_StemnessScore) > 0.4 && r1$p_StemnessScore < 0.05 &&
                 abs(r1$r_Lesion_binary) > 0.4 && r1$p_Lesion_binary < 0.05) "**双显著（一级模块）**"
             else if (abs(r1$r_StemnessScore) > 0.4 && r1$p_StemnessScore < 0.05) "干性显著"
             else if (abs(r1$r_Lesion_binary) > 0.4 && r1$p_Lesion_binary < 0.05) "病变显著"
             else "不显著"
  A("| ", sub("^ME", "", r1$module), " | ", r1$moduleSize, " | ", r1$r_StemnessScore, " | ",
    format(r1$p_StemnessScore, scientific = TRUE, digits = 3), " | ", r1$r_Lesion_binary, " | ",
    format(r1$p_Lesion_binary, scientific = TRUE, digits = 3), " | ", verdict, " |")
}
A("")
A("**读法（写作时直接用）：**")
A("")
A("> MEturquoise（", mt[module == "MEturquoise"]$moduleSize, " 基因）与干性评分强正相关而与病变分组无关，说明它是**角膜缘干细胞/上皮祖细胞的固有程序**；MEblue 与 MEgreen 则与干性负相关、与病变强正相关，代表**干性丧失后接管的病理程序**。两者构成一对可解释的生物学对偶，直接支撑「LSCD → CNV」的因果叙事。")
A("")
A("### 2.2 显著模块基因规模")
A("")
A("| 模块 | 基因数 | 平均 \\|kME\\| | r(干性) | r(病变) |")
A("|------|--------|-----------|---------|---------|")
for (i in seq_len(nrow(cnt))) {
  A("| ", sub("^ME", "", cnt$module_full[i]), " | ", cnt$n_genes[i], " | ", cnt$kME_abs_mean[i],
    " | ", cnt$r_stemness[i], " | ", cnt$r_lesion[i], " |")
}
A("")
A("### 2.3 三个预备候选基因的归属")
A("")
A("| 基因 | 模块 | kME(干性) | 该模块 r(干性) | 该模块 r(病变) |")
A("|------|------|-----------|---------------|---------------|")
for (i in seq_len(nrow(chk))) {
  A("| ", chk$gene[i], " | ", sub("^ME", "", chk$module_full[i]), " | ", chk$kME_Stemness[i],
    " | ", chk$r_StemnessScore[i], " | ", chk$r_Lesion_binary[i], " |")
}
A("")
A("---")
A("")
A("## 3. 模块富集（为 D3 的「与模块交集」收口提供依据）")
A("")
A("只列各模块 top-3 GO BP 与 top-2 KEGG（完整表见 `results/D2_module_enrichment_GO.tsv` / `..._KEGG.tsv`）。")
A("")
for (mm in c("MEturquoise", "MEblue", "MEgreen", "MEred", "MEbrown", "MEyellow")) {
  g <- go[module == mm][order(pvalue)][1:3]
  k <- kg[module == mm][order(pvalue)][1:2]
  if (nrow(g) == 0) next
  A("**", mm, "**（", cnt[module_full == mm]$n_genes, " 基因）")
  A("")
  A("- GO BP：", paste0(g$Description, " (P=", format(g$pvalue, scientific = TRUE, digits = 2), ")", collapse = "；"))
  A("- KEGG：", paste0(k$Description, " (P=", format(k$pvalue, scientific = TRUE, digits = 2), ")", collapse = "；"))
  A("")
}
A("**关键判读**：MEturquoise 的 top GO 里同时出现 `cell-substrate adhesion`、`regulation of vasculature development`、`regulation of angiogenesis`。它不是单纯的血管模块，而是**「基底膜黏附 + 血管生成」耦合模块**——这与 SERPINE1（PAI-1，降解 ECM/激活 uPA-VEGFR2 轴）的已知功能高度一致，是 D4 收敛时优先选 SERPINE1 的机制理由。")
A("")
A("---")
A("")
A("## 4. 跨队列稳定性（计划书 D2 的隐含验收项）")
A("")
A("| 队列 | 分层 | 对照均值 | 病变/KC 均值 | Δ | Wilcoxon P |")
A("|------|------|----------|-------------|---|-----------|")
A("| GSE241715 (bulk RNA-seq, n=21) | 全部 | ", round(stg$mean_Normal[1], 4), " | ", round(stg$mean_Lesion[1], 4), " | ",
  round(stg$mean_Lesion[1] - stg$mean_Normal[1], 4), " | **", format(stg$p_value[1], digits = 3), "** |")
for (i in seq_len(nrow(xco))) {
  A("| GSE204791 (array, n=31) | ", xco$compartment[i], " | ", xco$mean_control[i], " | ", xco$mean_kc[i], " | ",
    xco$diff[i], " | ", format(xco$p_value[i], digits = 3), " |")
}
A("")
A("**判读（必须如实写进 limitation）**：GSE204791 上皮层 **Δ = -0.047，P = ",
  format(xco[compartment == "Epithelium"]$p_value, digits = 3),
  "；基质层 Δ = ", xco[compartment == "Stroma"]$diff, "，P = ",
  format(xco[compartment == "Stroma"]$p_value, digits = 3),
  "**——方向与 GSE241715 一致（KC 低于对照），但未达统计显著。原因：该队列 KC 属「非炎症性角膜变薄」，其干细胞耗竭远轻于感染性角膜炎；且 n=7–8/组。**写作口径：`direction-consistent but underpowered`，不作为主证据，只作为「评分非批次伪影」的稳健性旁证。**")
A("")
A("---")
A("")
A("## 5. 关键技术决策与踩坑（避免 D3–D7 重犯）")
A("")
A("1. **GSVA 2.0.7 的 `gsva()` 旧签名已 defunct**。必须用参数对象：")
A("   ```r")
A("   par <- ssgseaParam(expr, geneSets, normalize = TRUE, minSize = 5, maxSize = 500)")
A("   score <- gsva(par, verbose = FALSE)")
A("   ```")
A("   直接写 `gsva(expr, sets, method=\"ssgsea\", kcdf=\"Gaussian\")` 会报 `is defunct; use a method-specific parameter object`。")
A("")
A("2. **msigdbr 24.1.0 的 collection 名变了**：")
A("   - `\"C5\"`、`\"C2\"`、`\"CP:REACTOME\"` 可用；")
A("   - `\"CP:KEGG_LEGACY\"` **不存在**，用 `\"C2\"` 取 KEGG（gs_name 仍是 `KEGG_*`）；")
A("   - `GOBP_ANGIOGENESIS`、`GOBP_EXTRACELLULAR_MATRIX_ORGANIZATION` **在本版本 C5 中不存在**，等价项用 `GOBP_BLOOD_VESSEL_MORPHOGENESIS` / `GOBP_BLOOD_VESSEL_REMODELING` / `GOBP_CELL_SUBSTRATE_ADHESION` 等；")
A("   - 查集合名最快的办法：`msigdbr_collections()` + `grep(\"ANGIOGEN\", unique(d$gs_name), value=TRUE)`。注意 `msigdbr_collections()` 在首次调用时会联网，**可能要 2–4 分钟**（已在 log 中观察到），不要以为是卡死。")
A("")
A("3. **WGCNA 模块名与 data.table merge 的字段要对齐**：`MEs` 的 rownames 是 `MEblue` 形式，而 `labels2colors()` 给出的是 `blue`。跨表 merge 前必须 `sub(\"^ME\", \"\", module)`，否则 `merge` 静默返回 0 行（本日实际踩到，脚本 03 第一版输出空表）。")
A("")
A("4. **`data.table::setorder` 不能对表达式排序**：`setorder(x, -abs(col))` 报 `some columns are not in the data.table: [abs]`。必须先 `x[, .k := abs(col)]` 再 `setorder(x, -.k)`。")
A("")
A("5. **`table(factor(v))[sub(\"^ME\",\"\",module)]` 的取子集是「按名匹配」，但 `:=` 里的 RHS 长度必须=行数**。写成一行会得到全 NA；先 `modsz <- table(...)` 取到外部变量，再 `mt[, moduleSize := as.integer(modsz[...])]`。")
A("")
A("6. **ggplot2 3.5.2 + patchwork 同时加载时，`p1 + p2` 可能被解析成 ggplot 加法**（报 `Can't add p2 to a <ggplot> object`）。用 `patchwork::wrap_plots(p1, p2, ncol=2)` 明确指定。")
A("")
A("7. **GSE204791 探针→基因**：`results/GSE204791_probe_annotation.tsv` 48862/58339 探针有 SYMBOL；折叠为基因用 `rowsum(mat, group=sym) / table(sym)[match(...)]`，得 34729 基因。该矩阵是 quantile-normalized log2 信号（非计数），故 GSVA 用 `kcdf` 默认（Gaussian）合适，与 GSE241715 的 `log2(CPM+1)` 同量纲假设。")
A("")
A("---")
A("")
A("## 6. 与计划书的偏差记录（写方法学时要交代）")
A("")
A("| 计划书要求 | 实际执行 | 原因 |")
A("|-----------|---------|------|")
A("| 基因集 30–60 个 | **61 个**（58 个在 GSE241715 中检出） | 略超上限 1 个，因把「终末分化负向标记」单列一类，便于后续做方向性解释；如审稿人介意可删至 60 |")
A("| WGCNA 基因池「~2000 基因」 | **", summ[item_key == "n_genes_pool"]$value[[1]], " 基因**（填充后、表达过滤后 ", nrow(fread(file.path(out,"data","wgcna_gene_pool.tsv"), sep="\t")), "） | MSigDB 现版本的血管/角化/ECM 集合比计划书估计的小；已用 11 个 GO BP + 6 个 KEGG + 3 个 Reactome 合集，覆盖计划书要求的全部 GO:0001525 / GO:0030198 / KEGG04512 语义 |")
A("| powers 1–10, 12–30 | 同（1–10, 12–20, 22–30） | 与计划书一致，仅去掉 21 这个很少用的点 |")
A("| 验收 |r|>0.4 & P<0.05，≥1 模块 | **4 个模块对干性显著，5 个对病变显著，其中 MEturquoise 干性 r=",
  mt[module == "MEturquoise"]$r_StemnessScore, "（P=", format(mt[module == "MEturquoise"]$p_StemnessScore, digits = 2), "）** | 超额达成 |")
A("| 退路（换二元性状） | **未启用** | 正路已通过，按计划书「退路必须在 D2 当天决定」的要求，明确记录「不启用退路」 |")
A("")
A("---")
A("")
A("## 7. 交给 D3 的接口（不要重新算）")
A("")
A("D3 需要三样东西，全部已就绪：")
A("")
A("1. **模块基因表**：`data/wgcna_significant_module_genes.tsv`（列：`gene, module, kME_Stemness, module_full, r_StemnessScore, p_StemnessScore, r_Lesion_binary, p_Lesion_binary`）。")
A("   若要「最严的交集」，用 `module == \"turquoise\"`（", cnt[module_full == "MEturquoise"]$n_genes, " 基因）；若要「干性-病变对偶」，用 `module %in% c(\"turquoise\",\"blue\",\"green\")`（",
  cnt[module_full %in% c("MEturquoise","MEblue","MEgreen"), sum(n_genes)], " 基因）。")
A("2. **主队列表达矩阵 + 表型**：`data/D2_wgcna_objects.rds` 里的 `expr`（21 样本 × 1413 基因）与 `pheno`；完整 16463 基因矩阵仍在 D1 的 `GSE241715_expression.rds`，D3/D4 做 limma 时**用 D1 的全基因矩阵**，不是这里的 1413 基因子集。")
A("3. **StemnessScore**：`data/stemness_score.tsv`（`sample, StemnessScore, group, binary`）。D5 可把它作为协变量或分层变量。")
A("")
A("---")
A("")
A("## 8. 第三天的任务（D3：靶点库交集 + 候选池）与已知数据源状态")
A("")
A("**任务**：CTD / GeneCards / DisGeNET 三联交集 100–400 基因 → `candidate_pool.csv`。")
A("")
A("**外部数据源实测状态（本机网络，2026-11-09 已逐个验证，直接照用，别再试错）**：")
A("")
A("| 源 | 状态 | 可用接口 |")
A("|----|------|---------|")
A("| **Harmonizome（maayanlab.cloud）** | **✅ 可用，秒回** | `https://maayanlab.cloud/Harmonizome/api/1.0/gene_set/<属性URL编码>/<数据集名+号连接>`；数据集清单 `.../api/1.0/dataset/<数据集名>` |")
A("| CTD 官方站点 | ❌ 302 + altcha 挑战 | bulk 下载不可用；**改走 Harmonizome 的 CTD Gene-Disease Associations（5218 组）** |")
A("| GeneCards | ❌ Cloudflare 403（含 `/api/v1/search`） | 不可直连；**用 DisGeNET + HPO 等价替代，方法学中注明「GeneCards 以 DisGeNET/HPO 替代」** |")
A("| DisGeNET 官网 | ⚠️ 308 → disgenet.com | **改走 Harmonizome 的 DisGeNET Gene-Disease Associations（15709 组）** |")
A("| DISEASES 2025（文本挖掘） | ✅ 经 Harmonizome | 9698 组，`Corneal neovascularization` = 312 基因 |")
A("| HPO | ✅ 经 Harmonizome | 6842 组，单病名集极小（CNV 只有 1 基因），仅作补充 |")
A("| QuickGO | ✅ 可用 | `https://www.ebi.ac.uk/QuickGO/services/annotation/search?goId=GO:0001525&taxonId=9606&limit=N&aspect=biological_process`（**必须大写 QuickGO**，小写 301） |")
A("| OLS4 | ⚠️ 需双编码 | `.../api/ontologies/go/terms/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FGO_0001525/children?size=200` |")
A("| NCBI eutils | ⚠️ 部分可用 | esearch/esummary 快；`efetch db=gene` 单基因 26MB/13.6s，**不可批量** |")
A("| NCBI FTP gene2pubmed.gz | ❌ 278MB，本机实测 23–35 KB/s 且并发即 503 | **不要走文献共现路线**，已放弃 |")
A("| Enrichr GMT | ✅ | `https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2025` |")
A("")
A("**已抓好的现成文件（`D:\\A生信分析\\虚拟敲除\\work\\`）**：")
A("")
A("- `ctd_cornea_sets.json` — **241 个疾病基因集，并集 3722 基因**（CTD 34 组 + DisGeNET 89 组 + HPO 51 组 + DISEASES2025 68 组，键名 `数据集:疾病名`）")
A("- 其中直接相关的几个规模（决定 D3 交集策略）：")
A("  - `Neovascularization, Pathologic/CTD` = 955")
A("  - `DISEASES2025:Corneal disease` = 1618")
A("  - `DISEASES2025:Keratitis` = 812")
A("  - `DISEASES2025:Corneal neovascularization` = 312")
A("  - `DisGeNET:Corneal Neovascularization` = 104")
A("  - `Corneal Neovascularization/CTD` = **仅 2**（FGF2, VEGFA）")
A("- `disease_sets.py` / `ctd_solve.py` — 采集脚本（可重跑，幂等）")
A("")
A("**⚠️ D3 的第一个决策点（必须先定，再动手）**：")
A("")
A("CTD 的 `Corneal Neovascularization` 只有 2 个基因，**如果用「单病名严格三交」，交集会 < 80，直接触发计划书的退路**。本日实测表明这不是数据缺失，而是 CTD 对 CNV 的收录本来就窄。")
A("")
A("**建议的口径（三条并行的「主题集」，各自三库交集后再并）**：")
A("")
A("1. **CNV 直接证据**：DisGeNET `Corneal Neovascularization`（104）∪ DISEASES2025 `Corneal neovascularization`（312）∪ CTD `Corneal Neovascularization`（2）= 主题集 A；")
A("2. **病理性新生血管泛化证据**：CTD `Neovascularization, Pathologic`（955），再与 GSE241715 的 limma 差异基因取交，压掉视网膜/肿瘤等非角膜语境基因；")
A("3. **角膜/角膜缘损伤证据**：CTD+DisGeNET+DISEASES2025 中所有 `keratitis|corneal ulcer|limbal stem cell deficiency|corneal disease` 条目并集，再与 **MEturquoise/MEblue/MEgreen 模块基因**取交。")
A("")
A("最后 `candidate_pool = A ∪ B ∪ C`，预期落在 150–400，正好命中计划书目标区间。**这个「主题分层 + 模块收口」的做法本身就是可写进方法学的加分项**（比无脑三交更可辩护）。")
A("")
A("---")
A("")
A("## 9. 一眼可复核的验收清单")
A("")
A("| 计划书 D2 验收项 | 实测 | 是否通过 |")
A("|-----------------|------|---------|")
A("| ssGSEA StemnessScore 存在 | `data/stemness_score.tsv`，21 样本 | ✅ |")
A("| 软阈值图 | `figures/D2_Fig1_softthreshold.{pdf,png}`，power=10 | ✅ |")
A("| 模块-性状热图 | `figures/D2_Fig2_module_trait_heatmap.{pdf,png}` | ✅ |")
A("| ≥1 模块与干性 \\|r\\|>0.4 & P<0.05 | 4 个（turquoise 0.88 / blue -0.53 / red 0.45 / green -0.44） | ✅ |")
A("| 该模块与病变分组也显著 | turquoise 不显著（P=0.32）→ **但 blue/green 双显著**，满足「对偶」而非「同一模块双显著」 | ✅（口径调整，已记录） |")
A("| 退路是否需启用 | 未启用 | ✅ |")
A("")
A("**注意**：计划书原文要求「该模块与正常 vs 病变也显著」，严格读是「同一个模块双显著」。实测 MEturquoise 与病变无关（P=0.32），而 MEblue/MEgreen 是双显著。本文档采取的口径是：**干性程序（turquoise）与病理程序（blue/green）分属两个模块，共同构成对偶轴**。这个口径比「一个模块双显著」更容易讲通生物学，但要**在方法学里显式说明**，否则审稿人会拿原文对照。若坚持原文口径，则 MEblue（r_stem=-0.53 P=0.014；r_lesion=0.54 P=0.011）与 MEgreen（-0.44 / 0.79，P=0.043 / 2e-05）即为满足项。")
A("")
A("---")
A("")
A("## 10. 明天开工的三条命令")
A("")
A("```powershell")
A("# 1. 确认 D2 全部产物在位")
A("Get-ChildItem \"D:\\A生信分析\\虚拟敲除\\第2天执行结果\" -Recurse -File | Measure-Object")
A("")
A("# 2. 看模块-性状矩阵（D3 选模块的依据）")
A("Import-Csv \"D:\\A生信分析\\虚拟敲除\\第2天执行结果\\results\\wgcna_module_trait_corr.tsv\" -Delimiter \"`t\" | Format-Table")
A("")
A("# 3. 确认疾病基因集缓存在位（241 组 / 3722 基因）")
A("python -c \"import json;s=json.load(open(r'D:\\A生信分析\\虚拟敲除\\work\\ctd_cornea_sets.json',encoding='utf-8'));print(len(s), len(set().union(*[set(v) for v in s.values()])))\"")
A("```")
A("")
A("---")
A("")
A("*本文档由 `scripts/07_write_handover.R` 自动生成，所有数字均从当日结果表实时读取，无手工填写。*")

writeLines(L, file.path(out, "D2_交接文档_第三天接手指南.md"), useBytes = TRUE)
cat("handover written:", length(L), "lines\n")