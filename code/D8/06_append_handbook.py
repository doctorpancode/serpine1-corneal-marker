# -*- coding: utf-8 -*-
"""D8-06 把 D8 执行内容追加进 D7 交接文档，并回写第 10 节的"后续可做"状态。"""
import io, os, re

DOC = r"D:/A生信分析/虚拟敲除/第7天执行结果/D7_交接文档_项目收尾.md"
OUT = r"D:/A生信分析/虚拟敲除/第8天执行结果/metadata/D8_handbook_appendix.md"

D8 = r"""
---

# D8 执行结果（扩网络 · 第三敲除 · 可编辑图）

> 本章为 D8 追加，承接第 10 节三条"后续可做"，现已全部完成。所有数值可回溯到
> `第8天执行结果/results/` 与 `第8天执行结果/logs/`，图在 `第8天执行结果/figures/`。

## D8-0 一句话结论

把网络从 494 扩到 947 个基因后，**"敲除 → 角质细胞分化程序被扰动"这一现象仍然成立，但它的特异性显著下降**：
947 个随机敲除里绝大多数都能扰动同一批上皮基因（Krt12 在 430/947 个敲除的 top-10 命中里出现），
而"敲除效应强度"与"该基因在网络里的出度"相关 ρ = 0.917。因此 D7 的结论要改写成
**"SERPINE1/TMSB4X 处在一个被网络拓扑放大的上皮分化模块上，二者方向相反"**，
而不能写成"只有敲这两个基因才会扰动上皮"。

## D8-1 扩网络：950 基因输入 → 947 基因网络

| 项目 | D7（494 网） | D8（947 网） |
|---|---|---|
| 目标基因数 | 500 | 950（`checkMemory` 上限 4503，取 0.85 安全系数内的目标值） |
| scQC 后存活 | 494 基因 × 2006 细胞 | **947 基因 × 3587 细胞** |
| 输入细胞 | 5000（下采样） | 5000（下采样自 8000 细胞的 QC 存活的 8994 基因 × 6657 细胞） |
| 网络非零边 | 169 562 | **600 604**（密度 66.97%） |
| 张量分解 norm explained | 47.7% | **56.4%** |
| 建网耗时 | — | 0.2 min（10 个网络，`nCores=6`） |

- 基因选择优先级：12 个核心候选 → D3 候选池（人→鼠符号映射）→ `VariableFeatures` → 其余，同组内按
  `tier`（`T1_triple_DEG` 优先）+ `n_source_sets` + 检出率排序。**Cdh2 经 scQC 被剔除**（检出率不足），
  所以最终核心基因是 11 个而非 12 个。
- 出/入度（`results/D8_network_edge_summary.tsv`）：Tmsb4x 833/841、Gja1 810/821、P4hb 804/811、
  Prdx2 797/804、Vim 793/804、S100a13 782/791、Vdac1 749/754、Sod1 745/748、Serpine1 729/737、
  Capg 709/717、**Xiap 362/329（明显偏低）**。
- 脚本：`scripts/01_build_network.R`；日志：`logs/01_build_network.txt`。
- **踩坑（已在脚本内绕过）**：官方 `scTenifoldKnk::scQC()` 的 MT 过滤写成 `X[mtGenes, ]`，
  在 `dgCMatrix` 上会 drop 成向量并抛
  `Error in base::colSums(x, na.rm = na.rm, dims = dims, ...) : 'x' must be an array of at least two dimensions`。
  D7 的 500 基因恰好不含 `MT-` 基因所以没触发，D8 的 950 个里有 1 个。脚本内改用逐条等价的
  `scQC_safe()`（`X[mt, , drop=FALSE]`）。

## D8-2 全基因敲除扫描（`transcriptomeWide` 路线）

- 脚本 `scripts/02_transcriptome_wide_scan.R`，直接调用官方
  `scTenifoldKnk::hkManifoldAlignment(WT, qcX, gKO = as.list(targets), t = 10)`，
  得到 `perturbationDistances` 947 KO × 947 genes 与 `perturbationDirections` 947 × 947。
- **耗时 0.03 min**（热核只需 `H` 的一行，批量敲除几乎免费），对比 D7 逐基因网重建是天壤之别。
- 输出：`results/D8_knockdown_spectrum.tsv`（长表 56.6 MB）、
  `results/D8_perturbation_distance_matrix.tsv`（18.5 MB）、`results/D8_ko_impact_rank.tsv`、
  `results/D8_top10_shared_targets.tsv`；原始对象 `data/D8_scan_raw.rds`。

**谱宽度 mean|Δ| 排名**（`D8_ko_impact_rank.tsv`，全 947 个敲除，核心基因的位次）：

| 基因 | mean\|Δ\| | 全库位次 /947 |
|---|---|---|
| S100a11 | 0.0712 | 1 |
| Hsp90ab1 | 0.0571 | 2 |
| mt-Co1 | 0.0519 | 3 |
| S100a6 | 0.0491 | 4 |
| **Tmsb4x** | **0.0450** | **5** |
| **Capg** | 0.00744 | 54 |
| **Prdx2** | 0.00532 | 73 |
| **Vdac1** | 0.00448 | 86 |
| **Gja1** | 0.00346 | 112 |
| **P4hb** | 0.00187 | 217 |
| **Vim** | 0.00153 | 270 |
| **S100a13** | 0.00118 | 337 |
| **Sod1** | 0.00065 | 484 |
| **Serpine1** | 0.00043 | 553 |
| **Xiap** | 6.3e-05 | 753 |

→ 只有 TMSB4X 是"全网前 5 的强扰动源"；**SERPINE1 落在第 553 位**，本身是个弱扰动基因。

**被扰动最多的是上皮/角化基因**（947 个敲除的 top-10 命中频次）：
Krt12 430、Lypd2 414、Mgarp 413、Slurp1 401、Aqp5 386、Ly6a 355、Lgals7 349、Aldh3a1 343、
Ly6e 333、Stmn1 322、Tmem176a 320、Smc4 307、Tubb5 302 …（`D8_top10_shared_targets.tsv`）

**与 D7 的一致性**（`D8_d7_d8_consistency.tsv`，D7 与 D8 网络共有 442 个基因）：

| 敲除 | ρ(mean\|Δ\| 排序) | top-10 重合基因 |
|---|---|---|
| Tmsb4x | 0.667 | Krt12, Aqp5, Gsta4, Upk3bl, Upk1b, Tmprss11d, Calml3（7 个） |
| Serpine1 | 0.412 | 无 |

→ TMSB4X 的靶点高度可复现；SERPINE1 的靶点在大网络里被背景淹没。

## D8-3 第三个（及第四、第五个）敲除

- 脚本 `scripts/03_third_ko.R`：在同一 947 基因网络上分别敲 Serpine1 / Tmsb4x / **P4hb / S100a13 / Vim**，
  每个 KO 走官方 `WT[g,] <- 0 → manifoldAlignment(d=2) → scTenifoldKnk::dRegulation(empiricalNull=TRUE, direction=...)`。
- 输出：`results/D8_third_ko.tsv`、`results/D8_third_ko_summary.tsv`、
  `results/D8_sig_overlap.tsv`、`results/D8_epithelial_program_score.tsv`。

**结果必须诚实表述**：**5 个敲除都只有 1 个显著基因，而且都是自身**；非目标显著基因数 = 0，
显著集两两 Jaccard 全为 0（`D8_sig_overlap.tsv`）。所以"第三敲除验证上皮程序"在**统计显著性层面没有得到支持**，
能支持的是**排序一致性**：

| 敲除 | top-5 非目标基因 | 上皮程序 18 基因平均 Z | 方向 |
|---|---|---|---|
| Tmsb4x | Krt5, S100a11, Krt12, Lypd2, Perp | **1.811** | 全部 down |
| P4hb | Krt5, Perp, Krt12, Lypd2, Calml3 | **1.728** | up |
| S100a13 | Krt5, Perp, Krt12, Lypd2, Calml3 | **1.583** | up |
| Serpine1 | Krt5, Perp, Krt14, Calml3, Upk3bl | **1.496** | up |
| Vim | Krt5, Perp, S100a11, Anxa2, Krt6a | **1.493** | up |

- 五个敲除的 top-5 都以 **Krt5 / Perp / Krt12** 开头 —— 这正是"上皮程序（角质化）被扰动"的保守性证据，
  但它是**排序层面**的，不是 p 值层面的。
- 方向仍可写：**只有 TMSB4X 敲除让上皮基因下调，其余（含 SERPINE1）都上调**。这与 D7 的"同一模块相反两端"一致，
  并且把它扩展到了 P4HB / S100A13 / VIM。

**中心性阴性对照（必须写进 Limitation）**：`results/D8_centrality_null_check.tsv` —
947 个敲除的"上皮基因平均 |Δ|"与该基因的网络出度 Spearman **ρ = 0.9169, p ≈ 0**。
即：网络里连接越多的基因，敲掉它越会"打到上皮"，与它的生物学身份无关。
`results/D8_epi_selectivity_ratio.tsv` 给出量化：上皮效应 / 其余基因效应（log2）的中位数，
背景基因 **1.150**（n = 936）vs 核心候选 **1.348**（n = 11）——两组几乎完全重叠。
核心基因里真正"偏上皮"的只有 Tmsb4x（2.74）、Prdx2（2.40）、P4hb（2.35）三个。

## D8-4 图（新增 4 张主图 + 1 张可编辑矢量图）

| 文件（`figures/main/`） | 内容 | 数据源 |
|---|---|---|
| `Fig10_knockout_impact_spectrum`（PDF/PNG/SVG） | 敲除影响谱热图：947 KO × 47 个有效读出基因，行内归一化 | `D8_scan_raw.rds` |
| `Fig11_impact_landscape`（PDF/PNG/SVG） | 扰动谱 PCA 景观，11 个核心候选高亮 | `D8_ko_impact_rank.tsv` |
| `Fig12_epithelial_vs_rest`（PDF/PNG/SVG） | 上皮选择性：log2(上皮效应 / 其余效应) 的分布与核心基因落点 | `D8_epi_selectivity_ratio.tsv` |
| `Fig13_third_knockout`（PDF/PNG/SVG） | P4hb / S100a13 / Vim 三个敲除的 top-12 差异调控 Z | `D8_third_ko.tsv` |

- **Fig10 画法说明（避免被问）**：每个 KO 的谱先按自身第 95 百分位缩放再截断到 1.5，
  这样比较的是"打在哪里"而不是"打得多重"；列只保留在 ≥10% 敲除里有信号的基因，
  否则会出现整片空白列（初版 947×120 的 z-score 版本几乎全白，已弃用）。
- **Fig11**：点 = 一个虚拟敲除，坐标来自其扰动谱的前两个主成分；颜色 = mean|Δ|。
  PC1 69.0% / PC2 10.1%。TMSB4X 单独落在 PC1 正端，SERPINE1 与 P4HB、S100a13、Vim 聚在 PC1 负端。
- **Fig12**：横轴为选择性 log2 比值，虚线为 0（上皮与其余等价）。核心基因并未整体位于背景分布右侧上沿，
  只有 Tmsb4x / Prdx2 / P4hb 明显偏离。**这张图是本章（也是全项目）最诚实的一张图**。
- **Fig13**：图注已明确写"经 empirical-null 校准后无非目标基因达显著，故展示排序"。
- `scripts/04_spectrum_figure.R`（Fig10–13 首版）+ `scripts/04b_fix_fig10_fig12.R`（Fig10/Fig12 定稿版）。

## D8-5 Fig9 机制图改为可编辑矢量版

- 产物：`figures/editable/Fig9_mechanism_model_editable.svg`（10 975 B）+ 同目录 PNG 预览
  + `Fig9_edit_guide.md`（图层与改字说明）。
- 结构：**9 个 Inkscape 图层**（`00 box backgrounds` / `01 source / injury` / `02 injury infiltrate` /
  `03 SERPINE1 / endothelium` / `04 TMSB4X / epithelium` / `05 vascular outcome` / `06 arrows & relations` /
  `07 quantitative evidence` / `08 legend`），**32 个 `<text>` 元素全部是活字**（字体
  `Arial, Helvetica, sans-serif`），10 条 `<path>` 箭头（`arrow` / `inhib` 两种 marker）。
- 画布 `760 × 440` 单位，物理尺寸 `200mm × 116mm`（单栏宽），白色底（`<rect fill="#FFFFFF">` 在
  `layer-boxes` 最底层）。
- 内容与 D7 版的**实质差别**：
  1. **改正逻辑**：D7 版把 SERPINE1 → TMSB4X 画成串联上下游，与"同一模块相反两端"的结论矛盾。
     新版把 TMSB4X 画成**与 SERPINE1 并行的上皮支**，两者之间用红色虚线 T 头标注
     "same module, opposite sign"；TMSB4X → 血管只保留红色虚线 "hypothesis only"。
  2. 加了 D8 证据行：`third knockout, D8 (947-gene network): P4hb / S100a13 / Vim KO → same keratinocyte
     module (Krt5, Perp, Krt12)`。
  3. 所有数字都带条件（AUCell 0.223 vs 0.004；↑393× neutrophil / ↑18× macrophage；ρ = 0.82 with
     macrophage；Krt12 Z 3.59 / 3.66），底部"Evidence"面板一行一句标注来源表。
- **渲染坑（必须记住）**：SVG 若**不带白色底矩形**，`rsvg::rsvg_png()` 渲染会留下深色底，
  黑色标题看上去像"消失"，排查了很久才定位。现在底白已写进脚本
  （`s = '    <rect ... fill="#FFFFFF"/>\n'.format(W, H)` 作为 `layer-boxes` 第一个元素）。
- 生成脚本：`scripts/05_fig9_editable.py`（纯字符串拼 SVG，可重复生成）；
  预览：`scripts/05b_render_fig9_png.R`（`rsvg::rsvg_png(..., width = 1520)`）。
  本机**没有 Inkscape / Illustrator**，`rsvg` 包通过代理安装成功
  （`scripts/05c_install_rsvg.R`，`HTTP_PROXY=http://127.0.0.1:7881`）。
- 投稿前仍要做的两件事（需要人手工）：
  1. 在 Illustrator / Inkscape 里打开 SVG，把字体统一成期刊要求（多数期刊要 Arial 或 Helvetica，
     `paint-order` 与 marker 需确认矢量展开正常）；
  2. 检查 `TMSB4X ⇢ vessel` 那条虚线在最终版式里是否要保留（它目前只是假说）。

## D8-6 新增 / 更新的 Limitation（承接第 8 节的八条，编号续 9–12）

9. **扩到 947 基因后，"扰动上皮程序"不再具有特异性**。947 个虚拟敲除中，Krt12 出现在 430 个敲除的
   top-10 命中里；上皮选择性比值的中位数在背景基因里已达 1.15。所以**不能写**"SERPINE1/TMSB4X 特异性地
   扰动上皮分化"，只能写"二者处于该模块的相反两端"。
10. **敲除效应强度与网络出度高度共线**（ρ = 0.9169, p ≈ 0）。`mean|Δ|` 排序本质上是"网络中心性排序"的
    代理，不能当功能重要性读。若要答审稿人，需要外部验证（如真实 KO 实验、siRNA 数据）而不是网络指标。
11. **D8 的 `dRegulation` 在 947 网里没有非目标显著基因**（5 个 KO 皆然）。显著基因少不是"数据干净"，
    而是 `empiricalNull` 校准把 Z 压得很低。Fig13 因此只展示排序，不做显著性声明。
12. **`transcriptomeWide` 路径的结果依赖 `t = 10` 的热核温度与网络本身**；`t` 越大越平滑。
    目前 `t` 用的是官方默认 `ma_heatT = 10` 与 `knockoutDirection(t = 5)`，两处不同是官方设定，
    跨 KO 比较时可以接受，但不能把 `perturbationDistance` 当有量纲的效应量。

## D8-7 复现方式（D8 部分）

```r
# 均在 R 4.4.2 下；Windows 上 Rscript 不在 PATH，用全路径
# C:\Program Files\R\R-4.4.2\bin\Rscript.exe，并设 R_MAX_VSIZE=12Gb（内存 15.2 GB，空闲约 7 GB）

setwd("D:/A生信分析/虚拟敲除/第8天执行结果/scripts")
source("01_build_network.R")        # 950 基因 → 947 基因网络，0.5 min，输出 data/D8_WT_network.rds
source("02_transcriptome_wide_scan.R")  # 947 个 KO 批量扫描，0.03 min
source("03_third_ko.R")             # 5 个 KO 逐基因 dRegulation，约 4 min
source("04_spectrum_figure.R")      # Fig10–Fig13（首版）
source("04b_fix_fig10_fig12.R")     # Fig10 / Fig12 定稿
execute("python", "05_fig9_editable.py")   # Fig9 SVG（Python 3.14: C:\Python314\python.exe）
source("05b_render_fig9_png.R")     # 用 rsvg 渲染 Fig9 PNG 预览
```

依赖增量：`rsvg`（本机新装）、`png`、`ggrepel`（已装）。其余与 D7 相同
（scTenifoldKnk 2.0.0 / scTenifoldNet 1.4.3 / ComplexHeatmap 2.22.0 / ggplot2 3.5.2 / svglite 2.2.2）。

## D8-8 目录

```
第8天执行结果/
├─ data/       D8_knk_input.rds, D8_knk_qcCounts.rds, D8_WT_network.rds, D8_scan_raw.rds
├─ results/    D8_network_gene_list.tsv, D8_network_edge_summary.tsv,
│              D8_knockdown_spectrum.tsv, D8_perturbation_distance_matrix.tsv,
│              D8_ko_impact_rank.tsv, D8_top10_shared_targets.tsv, D8_ko_top10_targets.tsv,
│              D8_third_ko.tsv, D8_third_ko_summary.tsv, D8_sig_overlap.tsv,
│              D8_epithelial_program_score.tsv, D8_epi_selectivity_ratio.tsv,
│              D8_centrality_null_check.tsv, D8_d7_d8_consistency.tsv
├─ figures/main/      Fig10–Fig13（各 PDF + PNG 300dpi + SVG）
├─ figures/editable/  Fig9_mechanism_model_editable.svg / .png, Fig9_edit_guide.md
├─ logs/       01_build_network.txt, 02_transcriptome_wide_scan.txt, 03_third_ko.txt,
│              04_spectrum_figure.txt, 04b_fix_fig10_fig12.txt
├─ metadata/   D8_network_meta.tsv, D8_handbook_appendix.md
└─ scripts/    01–05 全部脚本（含 05c_install_rsvg.R）
```

## D8-9 第 10 节三条的完成状态

| 原任务 | 状态 | 产物 |
|---|---|---|
| 1. 扩网络到 800–1000 + transcriptomeWide 批量扫描 + 影响谱热图 | **已完成** | 947 基因网络、947 KO 扫描、Fig10/Fig11/Fig12 |
| 2. 加第三个敲除（P4HB / S100A13）验证上皮程序保守性 | **已完成（结论需降级表述）** | Fig13、`D8_third_ko*.tsv`、Limitation 9–11 |
| 3. Fig9 机制图矢量可编辑版 | **已完成** | `figures/editable/Fig9_mechanism_model_editable.svg` + 编辑指南 |

D7 第 10 节剩余两条（GSE204791 人角膜单细胞验证、D2–D7 一页纸总汇总表）**未在 D8 处理**，仍需后续完成。
"""

def main():
    with io.open(DOC, "r", encoding="utf-8") as fh:
        doc = fh.read()
    if "\n# D8 执行结果" in doc:
        head, _, _ = doc.partition("\n# D8 执行结果")
        doc = head
        print("existing D8 section stripped before re-append")

    # 把第 10 节三条任务改成"已完成，见 D8"，保留其余待办
    old_10 = """## 10. 明日/后续可做

1. **扩网络**：把基因数从 494 提到 800–1000（官方上限），或改用 `transcriptomeWide = TRUE` + `hkManifoldAlignment`（热核批量对齐，可一次读多个敲除，内存更省），做全 500 基因的敲除扫描，画"敲除影响谱"热图。
2. **加第三个敲除**（P4HB 或 S100A13）验证"上皮程序被扰动"是否保守。
3. **把 Fig9 机制图做成矢量可编辑版**（Illustrator / Inkscape），投稿前替换。
4. **补 GSE204791 人角膜的单细胞验证**（若可用），把"巨噬驱动、内皮表达"在人类数据里再验一次。
5. **做一个 D2–D7 总汇总表**（一页纸：7 天 × 每天 1 张代表图 + 1 句结论），用于答辩首页。"""
    new_10 = """## 10. 明日/后续可做

1. ~~**扩网络**：把基因数从 494 提到 800–1000（官方上限），或改用 `transcriptomeWide = TRUE` + `hkManifoldAlignment`（热核批量对齐，可一次读多个敲除，内存更省），做全 500 基因的敲除扫描，画"敲除影响谱"热图。~~
   → **D8 已完成**（947 基因网络、947 个敲除的批量扫描、Fig10–Fig12）。**结论需降级**：效应与网络出度 ρ = 0.917，见 D8-6 与 Limitation 9–11。
2. ~~**加第三个敲除**（P4HB 或 S100A13）验证"上皮程序被扰动"是否保守。~~
   → **D8 已完成**（P4hb / S100a13 / Vim 三个）。统计显著性层面未获支持，仅排序层面保守（Krt5 / Perp / Krt12 稳定居前），见 D8-3、Fig13。
3. ~~**把 Fig9 机制图做成矢量可编辑版**（Illustrator / Inkscape），投稿前替换。~~
   → **D8 已完成**，`第8天执行结果/figures/editable/Fig9_mechanism_model_editable.svg`，9 个图层、32 个活字 `<text>`，见 D8-5。**仍待人工**：在本机外用 Illustrator 统一字体并确认 marker 矢量。
4. **补 GSE204791 人角膜的单细胞验证**（若可用），把"巨噬驱动、内皮表达"在人类数据里再验一次。（**未做**）
5. **做一个 D2–D7 总汇总表**（一页纸：7 天 × 每天 1 张代表图 + 1 句结论），用于答辩首页。（**未做**）
6. **新增（D8 引出）**：用真实敲除/干预数据（公共 siRNA/CRISPR 筛选）替代网络中心性指标，检验 SERPINE1/TMSB4X
   对上皮程序的效应是否独立于网络拓扑。这是当前最薄弱、也最可能被审稿人攻击的一环。"""
    if old_10 in doc:
        doc = doc.replace(old_10, new_10)
        print("section 10 updated")
    else:
        print("WARN: section 10 block not matched verbatim; left untouched")

    doc = doc.rstrip() + "\n" + D8
    with io.open(DOC, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(doc)
    with io.open(OUT, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(D8)
    print("doc bytes:", os.path.getsize(DOC))
    print("appendix bytes:", os.path.getsize(OUT))

if __name__ == "__main__":
    main()