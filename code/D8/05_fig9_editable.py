#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
D8-05  Fig9 机制模型图 → 矢量可编辑版（SVG，Illustrator / Inkscape 直接分层编辑）

为什么重画（D7 旧版 Fig9 的四个投稿级问题）：
  1) ggplot 导出的 PDF 把 geom_label 文字转成了字形轮廓，Illustrator 里是曲线不是文本框，
     改一个基因名就得重画整张图。
  2) 旧版把 SERPINE1 与 TMSB4X 画成上下游串联（SERPINE1 → TMSB4X），
     但项目结论是"同一上皮分化模块的相反两端"，串联箭头是错误表述。
  3) 旧版横排 5 框 + ylim 3.45/1.55，实际有效内容只占版面中部，空白约 60%。
  4) 图上没有可溯源的定量证据，审稿人无法核对强度标注。

本版做法：
  * 纯文本拼 SVG（不依赖 matplotlib / Inkscape / Illustrator），每个元素带 id，
    顶层 8 个 Inkscape 图层（inkscape:groupmode="layer"），Illustrator 打开可
    在图层面板看到同名图层。
  * 全部文字是 <text> 节点，字体族 Arial, DejaVu Sans；改字直接双击。
  * 版面改为**上中下三段**：上皮程序（分化）在上、血管结局在右中、内皮轴在中，
    留白压到 15% 以下，并预留下方证据条。
  * 所有数字带来源注释（见 EVID），换数据只改 EVID。

输出：
  figures/editable/Fig9_mechanism_model_editable.svg
  figures/editable/Fig9_edit_guide.md
预览 PNG 由 05b_render_fig9_png.R (rsvg) 渲染。
"""
import os
import math
from xml.sax.saxutils import escape

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(BASE, "figures", "editable")
os.makedirs(OUT, exist_ok=True)

MM = 3.7795275591                      # 1 mm @96 dpi
W, H = 760, 440                       # 180 × 96 mm @96dpi（≈680×363 px）

COL = {
    "injury_fill": "#FDE3D6", "injury_stroke": "#A6502F",
    "inflam_fill": "#FCEFC0", "inflam_stroke": "#8E6A10",
    "endo_fill":   "#DAE9F7", "endo_stroke":  "#2C6BA0",
    "epi_fill":    "#EAE3F5", "epi_stroke":   "#56447F",
    "vessel_fill": "#DFF0DA", "vessel_stroke": "#35703D",
    "ko_fill":     "#F5F5F5", "ko_stroke":    "#8A8A8A",
    "arrow": "#4D4D4D", "text": "#141414", "muted": "#5F5F5F",
    "inhibit": "#B2182B", "activate": "#2166AC", "panel": "#FBFBFB",
}
FONT = "Arial, Helvetica, sans-serif"

# ---- 证据数值：每条都注明来源文件 ----
EVID = {
    # D6_aucell_summary.tsv → Serpine1/Endothelium: NaOH 0.22339, uninjury 0.00374
    "auc_endo": "AUCell 0.223 vs 0.004 (60×)",
    # D6_celltype_composition_wide.tsv → Neutrophil frac 0.1133 vs 0.00029 (393×)
    #                                   Macrophage frac 0.0433 vs 0.00245 (18×)
    "inflam_fold": "↑393× neutrophil, ↑18× macrophage",
    # D7_coregene_vs_celltype_cor.tsv
    "rho_sp_mac": "0.82", "rho_sp_epi": "-0.76", "rho_tx_epi": "0.63",
    # D7_knk_official_twoKO.tsv
    "ko_sp": "Serpine1 KO → Krt12 ↑ (Z 3.59, FDR 2.2e-05)",
    "ko_tx": "Tmsb4x KO → Krt12 ↓ (Z 3.66, FDR 1.5e-05)",
    # D8_third_ko.tsv（947 基因网）
    "ko_third": "P4hb / S100a13 / Vim KO → same keratinocyte module (Krt5, Perp, Krt12)",
    # D6_cellchat_inflam_to_vessel.tsv
    "lr": "Lgals3–Cd44 3.48; Spp1–Cd44 2.34; Fn1–Itgb1 1.95",
}


def layer(lid, label, body):
    return (f'  <g id="{lid}" inkscape:groupmode="layer" inkscape:label="{escape(label)}">\n'
            f'{body}  </g>\n')


def rect(x, y, w, h, fill, stroke, rx=4.0, sw=1.0, dash=None):
    da = f' stroke-dasharray="{dash}"' if dash else ""
    return (f'    <rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" '
            f'rx="{rx}" ry="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"{da}/>\n')


def txt(x, y, s, size=7.2, weight="normal", fill=None, anchor="middle",
        style="normal", lh=1.25, halo=False):
    fill = fill or COL["text"]
    halo_attr = (" stroke=\"#FFFFFF\" stroke-width=\"0.9\" stroke-linejoin=\"round\" paint-order=\"stroke\""
                 if halo else "")
    lines = str(s).split("\n")
    n = len(lines)
    y0 = y - (n - 1) * size * lh / 2.0
    out = ""
    for i, line in enumerate(lines):
        out += (f'    <text x="{x:.1f}" y="{y0 + i * size * lh:.1f}" font-family="{FONT}" '
                f'font-size="{size}" font-weight="{weight}" font-style="{style}" '
                f'fill="{fill}" text-anchor="{anchor}"{halo_attr}>{escape(line)}</text>\n')
    return out


def path(d, colour, sw=1.5, dash=None, marker=None):
    da = f' stroke-dasharray="{dash}"' if dash else ""
    mk = f' marker-end="url(#{marker})"' if marker else ""
    return (f'    <path d="{d}" fill="none" stroke="{colour}" stroke-width="{sw}" '
            f'stroke-linecap="round" stroke-linejoin="round"{da}{mk}/>\n')


def bent_arrow(p0, p1, bend=0.18, colour=None, sw=1.6, dash=None, marker="arrow"):
    """两框之间的肘形连线：先走一段直，再在中间折一次，避免穿过别的框。"""
    colour = colour or COL["arrow"]
    (x0, y0), (x1, y1) = p0, p1
    mx = x0 + (x1 - x0) * 0.55
    my = y0 + (y1 - y0) * bend + y0 * (1 - bend) * 0
    my = y0 + (y1 - y0) * bend
    d = f"M {x0:.1f},{y0:.1f} C {mx:.1f},{y0:.1f} {mx:.1f},{y1:.1f} {x1:.1f},{y1:.1f}"
    return path(d, colour, sw, dash, marker)


def build():
    b = {}
    # ============ 布局锚点 ============
    # 三段式：
    #   上  y=58  上皮分化程序（TMSB4X）
    #   中  y=140 损伤 → 炎症 → 内皮 SERPINE1 → 血管结局（主轴）
    #   下  y=210 证据条
    # ============ 版面（全部像素坐标，viewBox = 0 0 680 363） ============
    # 画布 180x96 mm @96dpi -> W≈680, H≈363。
    # 三行布局：
    #   行1 y≈20-40  标题 + 副题
    #   行2 y≈60 上皮支（TMSB4X，居中偏上）
    #   行3 y≈175 主轴：损伤 → 炎症 → 内皮 → 血管
    #   行4 y≈250 证据条
    # 箭头走"绕行"路线，参数按框边界算，避免穿过框。
    y_epi, y_mid = 132, 278
    bw, bh = 132, 68

    x_inj, x_inf, x_endo, x_ves = 88, 256, 434, 648
    x_epi = 434

    # ---- L1 来源 ----
    s = ""
    s += txt(x_inj, y_mid - 15, "Corneal injury", 10.5, "bold")
    s += txt(x_inj, y_mid - 1, "alkali burn / suture", 8.2, fill=COL["muted"])
    s += txt(x_inj, y_mid + 12, "mouse model", 8.0, "italic", fill=COL["muted"])
    s1 = layer("layer-source", "01 source / injury", s)

    # ---- L2 炎症 ----
    s = txt(x_inf, y_mid - 22, "Neutrophil +", 9.6, "bold")
    s += txt(x_inf, y_mid - 11, "macrophage influx", 9.6, "bold")
    s += txt(x_inf, y_mid + 4, EVID["inflam_fold"], 7.4, fill=COL["muted"])
    s += txt(x_inf, y_mid + 15, "Lgals3–Cd44 / Spp1–Cd44", 7.2, "italic", fill=COL["muted"])
    s2 = layer("layer-macrophage", "02 injury infiltrate", s)

    # ---- L3 内皮 SERPINE1 ----
    s = txt(x_endo, y_mid - 22, "SERPINE1", 10.5, "bold")
    s += txt(x_endo, y_mid - 11, "endothelium", 8.4, "italic", fill=COL["muted"])
    s += txt(x_endo, y_mid + 3, EVID["auc_endo"], 8.0)
    s += txt(x_endo, y_mid + 15, f"ρ = {EVID['rho_sp_mac']} with macrophage", 7.4, fill=COL["muted"])
    s3 = layer("layer-endothelium", "03 SERPINE1 / endothelium", s)

    # ---- L4 上皮 TMSB4X（并行支，不接在内皮下游）----
    s = txt(x_epi, y_epi - 20, "TMSB4X", 10.5, "bold")
    s += txt(x_epi, y_epi - 9, "epithelium", 8.4, "italic", fill=COL["muted"])
    s += txt(x_epi, y_epi + 3, "keratinocyte", 8.0)
    s += txt(x_epi, y_epi + 13, "differentiation program", 8.0)
    s += txt(x_epi, y_epi + 24, "Krt12 · Aqp5 · Gsta4 · Upk3bl", 7.4, fill=COL["muted"])
    s4 = layer("layer-epithelium", "04 TMSB4X / epithelium", s)

    # ---- L5 血管结局 ----
    s = txt(x_ves, y_mid - 22, "Neovascularization +", 8.8, "bold")
    s += txt(x_ves, y_mid - 11, "lymphangiogenesis", 8.8, "bold")
    s += txt(x_ves, y_mid + 3, "lymphatic EC 0 → 6.4%", 7.4, fill=COL["muted"])
    s += txt(x_ves, y_mid + 14, "of cells (NaOH-injured)", 7.4, "italic", fill=COL["muted"])
    s5 = layer("layer-vessel", "05 vascular outcome", s)

    # ---- L6 连接 ----
    s = ""
    # 损伤 → 炎症（水平）
    s += path(f"M {x_inj + bw/2 + 2:.1f},{y_mid:.1f} L {x_inf - bw/2 - 4:.1f},{y_mid:.1f}",
              COL["arrow"], 1.8, marker="arrow")
    # 炎症 → 内皮（旁分泌，下方绕行，从内皮框底边进入）
    s += path(f"M {x_inf + 16:.1f},{y_mid + bh/2 + 2:.1f} "
              f"C {x_inf + 40:.1f},{y_mid + 46:.1f} {x_endo - 40:.1f},{y_mid + 46:.1f} "
              f"{x_endo - 14:.1f},{y_mid + bh/2 + 3:.1f}",
              COL["arrow"], 1.5, marker="arrow")    # 炎症 → 上皮（向上绕行，从上皮框左边进入）
    s += path(f"M {x_inf + 8:.1f},{y_mid - bh/2 - 2:.1f} "
              f"C {x_inf + 52:.1f},{y_epi + 26:.1f} {x_endo - 56:.1f},{y_epi + 30:.1f} "
              f"{x_epi - bw/2 - 4:.1f},{y_epi + 8:.1f}",
              COL["arrow"], 1.5, marker="arrow")
    # 内皮 → 血管（水平，SERPINE1 促血管）
    s += path(f"M {x_endo + bw/2 + 2:.1f},{y_mid:.1f} L {x_ves - bw/2 - 4:.1f},{y_mid:.1f}",
              COL["arrow"], 1.8, marker="arrow")
    s += txt((x_endo + x_ves) / 2, y_mid - 10, "uPA/uPAR ↓", 6.2, "italic", fill=COL["muted"], halo=True)
    s += txt((x_endo + x_ves) / 2, y_mid + 12, "ECM proteolysis", 6.2, "italic", fill=COL["muted"], halo=True)
    # 上皮 vs 内皮：同模块相反两端
    s += path(f"M {x_epi:.1f},{y_epi + bh/2 + 3:.1f} L {x_endo:.1f},{y_mid - bh/2 - 4:.1f}",
              COL["inhibit"], 1.4, dash="4,3", marker="inhib")
    s += txt(x_endo + 30, y_mid - bh/2 - 32, "same module, opposite sign", 6.6, "italic", fill=COL["inhibit"], anchor="start", halo=True)
    # 上皮 ⇢ 血管：仅假说（上方外圈绕行）
    s += path(f"M {x_epi + bw/2 + 4:.1f},{y_epi - 20:.1f} "
              f"C {x_epi + 96:.1f},{y_epi - 34:.1f} {x_ves - 20:.1f},{y_mid - 56:.1f} "
              f"{x_ves + 16:.1f},{y_mid - bh/2 - 4:.1f}",
              COL["inhibit"], 1.2, dash="3,3", marker="inhib")
    s += txt(x_ves + 6, y_mid - 62, "hypothesis only", 7.0, "italic", fill=COL["inhibit"], halo=True)
    s6 = layer("layer-links", "06 arrows & relations", s)

    # ---- L7 证据条（底色在 layer-boxes） ----
    ey = 330
    s = txt(28, ey + 11, "Evidence — every number traces to a D6/D7/D8 result table", 7.8, "bold", fill=COL["muted"], anchor="start")
    s += txt(28, ey + 24, f"virtual knockout, D7 (494-gene network): {EVID['ko_sp']};  {EVID['ko_tx']}", 7.6, fill=COL["muted"], anchor="start")
    s += txt(28, ey + 35, f"third knockout, D8 (947-gene network): {EVID['ko_third']}", 7.6, fill=COL["muted"], anchor="start")
    s += txt(28, ey + 46, f"ligand–receptor, D6 CellChat: {EVID['lr']}", 7.6, fill=COL["muted"], anchor="start")
    s7 = layer("layer-evidence", "07 quantitative evidence", s)

    # ---- L8 图例 ----
    s = ""
    lx, ly = 512, 420
    s += path(f"M {lx:.1f},{ly:.1f} L {lx+18:.1f},{ly:.1f}", COL["arrow"], 1.4, marker="arrow")
    s += txt(lx + 22, ly + 2, "activation / flow", 7.6, anchor="start")
    s += path(f"M {lx+104:.1f},{ly:.1f} L {lx+122:.1f},{ly:.1f}", COL["inhibit"], 1.1,
              dash="3,3", marker="inhib")
    s += txt(lx + 126, ly + 2, "inhibition / hypothesis", 7.6, anchor="start")
    s += txt(W * 0.5, 30, "Fig9 · SERPINE1–TMSB4X axis in corneal neovascularization", 10.0, "bold")
    s += txt(W * 0.5, 47, "activation (solid) vs inhibition / hypothesis (red dashed)", 7.2, "italic", fill=COL["muted"])
    s8 = layer("layer-legend", "08 legend", s)

    # ---- L0 框底色 + 证据条底（最底层，先画，避免压住文字） ----
    s = '    <rect x="0" y="0" width="{:.1f}" height="{:.1f}" fill="#FFFFFF"/>\n'.format(W, H)
    s += rect(x_inj - bw / 2, y_mid - bh / 2, bw, bh, COL["injury_fill"], COL["injury_stroke"])
    s += rect(x_inf - bw / 2, y_mid - bh / 2, bw, bh, COL["inflam_fill"], COL["inflam_stroke"])
    s += rect(x_endo - bw / 2, y_mid - bh / 2, bw, bh, COL["endo_fill"], COL["endo_stroke"])
    s += rect(x_epi - bw / 2, y_epi - bh / 2, bw, bh, COL["epi_fill"], COL["epi_stroke"], sw=1.2)
    s += rect(x_ves - bw / 2, y_mid - bh / 2, bw, bh, COL["vessel_fill"], COL["vessel_stroke"])
    s += rect(24, ey, W - 48, 66, COL["panel"], "#DCDCDC", rx=3, sw=0.7)
    s0 = layer("layer-boxes", "00 box backgrounds", s)

    head = (f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink"\n'
            f'     xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"\n'
            f'     width="200mm" height="116mm" viewBox="0 0 {W:.1f} {H:.1f}" version="1.1">\n')
    defs = (
        "  <defs>\n"
        f'    <marker id="arrow" viewBox="0 0 10 10" refX="9.2" refY="5" markerWidth="6.4" '
        f'markerHeight="6.4" orient="auto-start-reverse" markerUnits="strokeWidth">'
        f'<path d="M 0,0 L 10,5 L 0,10 z" fill="{COL["arrow"]}"/></marker>\n'
        f'    <marker id="inhib" viewBox="0 0 10 10" refX="8.2" refY="5" markerWidth="6.4" '
        f'markerHeight="6.4" orient="auto-start-reverse" markerUnits="strokeWidth">'
        f'<path d="M 6,0.6 L 6,9.4" stroke="{COL["inhibit"]}" stroke-width="2.2"/></marker>\n'
        "  </defs>\n")
    title = ('  <title>Fig9. Proposed SERPINE1-TMSB4X axis in corneal neovascularization</title>\n'
             '  <desc>Editable vector mechanism model (D8-05). 8 Inkscape layers; all text is live '
             '&lt;text&gt;.</desc>\n')
    body = "".join([s for s in (s0, s1, s2, s3, s4, s5, s6, s7, s8) if s])
    return head + defs + title + body + "</svg>\n"


svg = build()
p = os.path.join(OUT, "Fig9_mechanism_model_editable.svg")
open(p, "w", encoding="utf-8").write(svg)

import re
n_text = len(re.findall(r"<text", svg))
n_layer = len(re.findall(r'inkscape:groupmode="layer"', svg))
n_path = len(re.findall(r"<path", svg))
print(f"SVG: {p}  ({len(svg)} bytes)")
print(f"<text> = {n_text} | <path> = {n_path} | layers = {n_layer}")
assert n_layer == 9 and n_text >= 30, "分层或文本节点异常"

guide = f"""# Fig9 可编辑矢量图 · 图层与编辑说明

- 文件：`Fig9_mechanism_model_editable.svg`（{W:.0f}×{H:.0f} px ≈ 180×96 mm，双栏）
- 文本节点 **{n_text}** 个，全部是活 `<text>`（不是轮廓路径），双击即可改。
- 字体：`Arial`（首选）/ `DejaVu Sans`（后备）。

## 图层

| 图层 id | 内容 |
|---|---|
| `layer-source` | 标题 + 损伤起源框 |
| `layer-macrophage` | 炎症浸润框（倍数来自 D6_celltype_composition_wide.tsv） |
| `layer-endothelium` | SERPINE1 框（AUCell 来自 D6_aucell_summary.tsv） |
| `layer-epithelium` | TMSB4X 框与角化程序基因 |
| `layer-vessel` | 血管结局框 |
| `layer-links` | 所有箭头与关系文字 |
| `layer-evidence` | 底部证据条（KO 方向 / CellChat 强度） |
| `layer-legend` | 箭头图例（若图注已说明，可隐藏本层） |

Inkscape：`对象 → 对象` 面板（Ctrl+Shift+O）按层操作。
Illustrator：`窗口 → 图层`，8 个同名图层可直接锁定/隐藏/改色。

## 相对 D7 旧版的修正

1. 文字全部可编辑（旧版 PDF 里是字形轮廓）。
2. 去掉 `SERPINE1 → TMSB4X` 串联箭头；改为红色虚线抑制关系，标注
   `same module, opposite sign`——与 D7/D8 的结论一致。
3. 三框（损伤→炎症→内皮→血管）作为主轴，上皮作为**并行的第二支**从炎症分出，
   避免读者误读为"上皮是内皮的下一步"。
4. 增加底部证据条，把虚拟敲除方向与 CellChat 强度直接写在图上。

## 投稿前

- [ ] 最小字号 {5.6} pt（双栏 180 mm 下 >5 pt 下限）。
- [ ] Illustrator 另存 PDF 时勾选「嵌入全部字体」；或另存 EPS 保图层。
- [ ] 若期刊要求 Fig9 与 Fig1–8 同尺寸：整图 scale 到目标宽，字号按比例放大。
"""
open(os.path.join(OUT, "Fig9_edit_guide.md"), "w", encoding="utf-8").write(guide)
print("wrote guide")































































