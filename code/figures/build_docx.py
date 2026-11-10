# -*- coding: utf-8 -*-
"""把 figs_v2/out 的图件与表1 插入 草稿1.0.docx → 草稿1.1_带图.docx

要点：
- 先用「原始段落对象」捕获锚点，再做插入，避免插入后索引漂移（上一版 表1 曾插错位置）。
- 图：锚点 → 图片段 → 粗体图题 → 小字图注。
- 表1：锚点 → 表题（与表保持同页）→ 三线表（分节行跨列合并）→ 表注。
"""
import os
import csv
from docx import Document
from docx.shared import Pt, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

BASE = r"D:\A生信分析\虚拟敲除\文章"
OUT = os.path.join(BASE, "figs_v2", "out")
SRC = os.path.join(BASE, "草稿1.0.docx")
DST = os.path.join(BASE, "草稿1.1_带图.docx")

CAPS = {
    "Fig1_design": (
        "图1 研究设计与分析流程",
        "(a) 六阶段分析流程及各阶段的输入与输出；(b) WGCNA 6 个模块与干性评分、病损状态的 "
        "Pearson 相关（*P<0.05，**P<0.01，***P<0.001，****P<0.0001，ns 为不显著）；"
        "(c) 三个候选基因池的交集与并集基因数。"),
    "Fig2_bulk": (
        "图2 bulk 队列的表达梯度与差异表达",
        "(a) GSE241715 四组 SERPINE1 表达（均值±标准差，n=7/4/5/5），呈 5.4 倍梯度，"
        "Kruskal–Wallis P=9.4×10⁻⁴；(b) 10 个核心基因的组间均值热图（色标为行 z-score，格内为原值）；"
        "(c) 病损与正常比较的差异表达（limma；虚线为 |log₂FC|=0.5 与 adjP=0.05，共 1 350 个差异基因）。"),
    "Fig3_singlecell": (
        "图3 单细胞图谱与损伤后的细胞组成重塑",
        "(a) 15 873 个细胞的 UMAP 与 11 类细胞注释；(b) 两组各细胞类型占比（标注为碱烧伤/未损伤，%）；"
        "(c) 占比的对数倍数。"),
    "Fig4_localization": (
        "图4 信号来源与表达位点的空间分离",
        "(a) SERPINE1 在各细胞类型的平均 AUCell 活性（标注为碱烧伤/未损伤）；"
        "(b) SERPINE1 与各细胞类型丰度的 Spearman 相关（21 例）；(c) TMSB4X 的同类相关；"
        "(d) NNLS 解卷积估计的四组细胞组成。"),
    "Fig5_communication": (
        "图5 免疫浸润、细胞通讯与通路富集",
        "(a) 10 个 adjP<0.05 的免疫细胞类型浸润分数差（*P<0.05）；(b) CellChat 通路强度倍数（标注比值；"
        "SPP1 与 TENASCIN 在未损伤侧为 0，比值仅为形式表达）；(c) 指向血管内皮的配体–受体对"
        "（Neu 中性粒、Mac 巨噬、Endo 内皮、LEC 淋巴管内皮）；(d) KEGG 与 Reactome 富集"
        "（气泡大小为基因数，颜色为 −log₁₀ adjP）。"),
    "Fig6_perturbation": (
        "图6 虚拟敲除的扰动幅度、网络位置与上皮选择性",
        "(a) Serpine1 与 Tmsb4x 敲除对 4 个上皮基因的 Z 值；(b) 10 个公共下游靶基因及其受扰动 KO 数；"
        "(c) 947 次敲除的出度与上皮 18 基因平均 |Δ|（双对数坐标，Spearman ρ=0.917）；"
        "(d) 剔除出度后的残差（括号内为残差排名）；(e) 上皮/非上皮效应比分布（候选 11 个 vs 背景 936 个）；"
        "(f) 出度三分位层内的相关。"),
    "Fig7_nuclear_effect": (
        "图7 核效应量级：SERPINE1 是标志而非节点",
        "(a) 18 个上皮分化基因在三个敲除下的 Δ；(b) 上皮程序平均 |Δ| 的量级差（480×）；"
        "(c) 逐基因 SERPINE1 Δ 与 TMSB4X Δ 的关系；(d) 947 次敲除的平均 |Δ| 分布与三个核心基因的位置。"),
    "S1_external_core": ("图S1 GSE58291 三个核心基因的三组表达",
                         "柱高为组均值、误差棒为标准差；标注为 BK 组相对 Control 的 log₂FC、"
                         "BH 校正 P 值与三组趋势 ρ。"),
    "S2_gradient": ("图S2 GSE58291 的病损程度梯度",
                    "(a) 三基因在对照、FK、BK 三组的均值折线（误差棒为 SD，图例内为趋势 ρ）；"
                    "(b) 相对对照组的 Δlog₂ 表达。SERPINE1 与 P4HB 单调上升，GJA1 下降后回升。"),
    "S3_five_ko": ("图S3 五个基因虚拟敲除的共享受扰靶点",
                   "(a) 靶点被几个敲除共享的分布（上图）与 10 个共享靶点的 Z 值（下图，格内为 Z，"
                   "右侧为该基因被共享的敲除数）；(b) 五个敲除对 18 个上皮基因的平均 |Z|。"),
    "S4_multiko": ("图S4 多基因联合敲除的上皮效应叠加",
                   "(a) 实测组合效应与单敲平均（菱形）的对比，右侧为实测/预测比值；"
                   "(b) 实测对线性预测的散点，虚线为 y=x，位于线上方表示超加性。"),
    "S5_spectrum_shift": ("图S5 归一化前后的影响谱重排",
                          "(a) 947 个敲除按上皮平均 |Δ| 的降序谱（红点为 11 个核心基因）；"
                          "(b) 行 L1 归一化前后的排名位移（空心为原始、实心为归一化）。"),
    "S6_multiko_genes": ("图S6 多敲除的逐基因效应矩阵",
                         "(a) 18 个上皮基因在三个敲除下的 Δ（色标为列内归一化，格内为原始 Δ）；"
                         "(b) 三敲绝对效应与符号一致性。"),
    "S7_GSE326016": ("图S7 GSE326016 的上皮程序得分",
                     "(a) 两组逐样本得分（横线为均值，误差棒为 SD，Welch t 检验 P=0.110）；"
                     "(b) 按得分排名配对的样本变化。"),
    "S8_rank_shift": ("图S8 行归一化前后的排名变化",
                      "(a) 11 个核心基因在原始与行归一化度量下的排名（蓝为前移、红为后退）；"
                      "(b) 三种度量与出度的 Spearman ρ。"),
}

ORDER = ["Fig1_design", "Fig2_bulk", "Fig3_singlecell", "Fig4_localization",
         "Fig5_communication", "Fig6_perturbation", "Fig7_nuclear_effect"]
SUPP = ["S1_external_core", "S2_gradient", "S3_five_ko", "S4_multiko",
        "S5_spectrum_shift", "S6_multiko_genes", "S7_GSE326016", "S8_rank_shift"]

ANCHOR = {"Fig1_design": 39, "Fig2_bulk": 41, "Fig3_singlecell": 40,
          "Fig4_localization": 45, "Fig5_communication": 48,
          "Fig6_perturbation": 56, "Fig7_nuclear_effect": 57}
T1_ANCHOR = 48

CN_FONT = "宋体"


def set_run_font(run, size=9, bold=False, font=CN_FONT):
    run.font.size = Pt(size)
    run.font.bold = bold
    run.font.name = font
    rpr = run._element.get_or_add_rPr()
    rf = rpr.find(qn("w:rFonts"))
    if rf is None:
        rf = OxmlElement("w:rFonts")
        rpr.append(rf)
    rf.set(qn("w:eastAsia"), font)


def new_par_after(par, doc):
    p = doc.add_paragraph()
    par._p.addnext(p._p)
    return p


def add_caption_after(par, doc, title, note):
    p = new_par_after(par, doc)
    p.alignment = WD_ALIGN_PARAGRAPH.LEFT
    pf = p.paragraph_format
    pf.space_before = Pt(6)
    pf.space_after = Pt(2)
    pf.line_spacing = 1.15
    pf.keep_with_next = True
    r = p.add_run(title)
    set_run_font(r, 9, bold=True)
    p2 = new_par_after(p, doc)
    p2.paragraph_format.space_before = Pt(0)
    p2.paragraph_format.space_after = Pt(14)
    p2.paragraph_format.line_spacing = 1.15
    r2 = p2.add_run("注：" + note if not note.startswith("(") else note)
    set_run_font(r2, 8.5, bold=False)
    return p, p2


def add_picture_after(par, doc, path, width_in=6.25):
    p = new_par_after(par, doc)
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_before = Pt(10)
    p.paragraph_format.space_after = Pt(2)
    p.paragraph_format.keep_with_next = True
    p.add_run().add_picture(path, width=Inches(width_in))
    return p


def _edges(tcPr, spec):
    borders = OxmlElement("w:tcBorders")
    for edge, sz in spec.items():
        e = OxmlElement("w:" + edge)
        e.set(qn("w:val"), "single" if sz else "none")
        e.set(qn("w:sz"), str(sz or 0))
        e.set(qn("w:space"), "0")
        e.set(qn("w:color"), "000000")
        borders.append(e)
    tcPr.append(borders)


def three_line_table(doc, rows, widths=(2.75, 1.05, 1.45, 0.75)):
    ncol = len(rows[0])
    t = doc.add_table(rows=len(rows), cols=ncol)
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    t.autofit = False
    tblPr = t._tbl.tblPr
    layout = OxmlElement("w:tblLayout")
    layout.set(qn("w:type"), "fixed")
    tblPr.append(layout)
    borders = OxmlElement("w:tblBorders")
    for edge in ("top", "left", "bottom", "right", "insideH", "insideV"):
        e = OxmlElement("w:" + edge)
        e.set(qn("w:val"), "none")
        e.set(qn("w:sz"), "0")
        borders.append(e)
    tblPr.append(borders)
    for j, w in enumerate(widths[:ncol]):
        t.columns[j].width = Inches(w)

    for i, row in enumerate(rows):
        is_group = (0 < i < len(rows) - 1 and row[0].strip()
                    and not any(c.strip() for c in row[1:]))
        t.rows[i].height = None
        trPr = t.rows[i]._tr.get_or_add_trPr()
        cant = OxmlElement("w:cantSplit")
        trPr.append(cant)
        if i == 0:  # 表头跨页重复
            hdr = OxmlElement("w:tblHeader")
            trPr.append(hdr)
        if is_group:
            merged = t.cell(i, 0).merge(t.cell(i, ncol - 1))
            cells = [merged]
        else:
            cells = [t.cell(i, j) for j in range(ncol)]
        for j, cell in enumerate(cells):
            if len(cells) > 1:
                cell.width = Inches(widths[j])
            cell.text = ""
            p = cell.paragraphs[0]
            p.alignment = (WD_ALIGN_PARAGRAPH.CENTER if j and not is_group
                           else WD_ALIGN_PARAGRAPH.LEFT)
            p.paragraph_format.space_before = Pt(5 if is_group else 1.5)
            p.paragraph_format.space_after = Pt(1.5)
            text = row[0] if is_group else row[j]
            r = p.add_run(str(text))
            set_run_font(r, 8.5, bold=(i == 0 or is_group))
            tcPr = cell._tc.get_or_add_tcPr()
            spec = {}
            if i == 0:
                spec = {"top": 12, "bottom": 6}
            elif i == len(rows) - 1:
                spec = {"bottom": 12}
            if spec:
                _edges(tcPr, spec)
    return t


def replace_text(par, old, new):
    for r in par.runs:
        if old in r.text:
            r.text = r.text.replace(old, new)
            return True
    if old not in par.text:
        return False
    txt = par.text.replace(old, new)
    for r in par.runs:
        r.text = ""
    par.runs[0].text = txt
    return True


def main():
    doc = Document(SRC)
    ps = doc.paragraphs
    # 先捕获锚点对象（插入会改变索引）
    anchors = {k: ps[v] for k, v in ANCHOR.items()}
    a_t1 = ps[T1_ANCHOR]

    if replace_text(ps[41], "正常 7 例、病毒性 5 例、细菌真菌性 5 例、圆锥角膜 4 例",
                    "正常 7 例、圆锥角膜 4 例、病毒性 5 例、细菌/真菌性 5 例"):
        print("修订 [41] 组别顺序")
    if replace_text(ps[39], "联合分析（GSE275678", "联合分析（图1；GSE275678"):
        print("补 图1 引用")
    if replace_text(ps[35], "主图以 300 dpi 的 PDF 与 PNG 双格式导出，示意图另存可编辑 SVG。",
                    "主图以 400 dpi 的 PNG 与 PDF 双格式导出。"):
        print("修订 [35] 导出规格")

    # 图5 先插，表1 后插 → 文档顺序为 表1、图5
    def insert_fig(key, anchor):
        path = os.path.join(OUT, key + ".png")
        if not os.path.exists(path):
            print("!! 缺图", path)
            return anchor
        title, note = CAPS[key]
        _t, tail = add_caption_after(anchor, doc, title, note)
        add_picture_after(anchor, doc, path)
        print("插入", key)
        return tail

    insert_fig("Fig1_design", anchors["Fig1_design"])
    t2 = insert_fig("Fig2_bulk", anchors["Fig2_bulk"])
    insert_fig("Fig3_singlecell", t2)          # 图3 紧随图2，保证图号顺序
    insert_fig("Fig4_localization", anchors["Fig4_localization"])
    insert_fig("Fig5_communication", anchors["Fig5_communication"])
    insert_fig("Fig6_perturbation", anchors["Fig6_perturbation"])
    insert_fig("Fig7_nuclear_effect", anchors["Fig7_nuclear_effect"])

    tsv = os.path.join(OUT, "Table1_diagnostic.tsv")
    if os.path.exists(tsv):
        with open(tsv, encoding="utf-8-sig") as f:
            rows = [r for r in csv.reader(f, delimiter="\t")]
        p_title, p_note = add_caption_after(
            a_t1, doc, "表1 诊断模型性能与共线性诊断",
            "AUC 为交叉验证曲线下面积；P 为 Hosmer–Lemeshow 校准检验；OR 为单因素 logistic 回归；"
            "VIF 为方差膨胀因子；— 表示该项不适用。")
        tbl = three_line_table(doc, rows)
        p_title._p.addnext(tbl._tbl)
        p_note.paragraph_format.space_before = Pt(4)
        print("插入表1")

    last = doc.paragraphs[-1]
    h = new_par_after(last, doc)
    h.paragraph_format.space_before = Pt(18)
    h.paragraph_format.keep_with_next = True
    r = h.add_run("补充材料")
    set_run_font(r, 12, bold=True)
    tail = h
    for key in SUPP:
        path = os.path.join(OUT, key + ".png")
        if not os.path.exists(path):
            print("!! 缺图", path)
            continue
        title, note = CAPS[key]
        prev = tail
        _t, tail = add_caption_after(prev, doc, title, note)
        add_picture_after(prev, doc, path)   # 插在 prev 之后 → 落到图题之前
        print("追加", key)

    doc.save(DST)
    print("saved", DST, "段落数", len(doc.paragraphs))


if __name__ == "__main__":
    main()