# -*- coding: utf-8 -*-
"""表1  诊断模型性能、单因素 OR 与共线性诊断
输出：out/Table1_diagnostic.tsv（可直接复制进 Word）+ out/Table1_preview.png（三线表预览）
数字全部来自第 5 天结果表，不修改任何数值。
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import figstyle as S

OUT = S.OUT
auc = S.load(5, "D5_model_ROCAUC.tsv")
orr = S.load(5, "D5_univariable_OR.tsv").set_index("gene")
vif = S.load(5, "D5_vif.tsv")

MODEL_NAME = {"M1_TMSB4X": "M1 TMSB4X", "M2_SERPINE1": "M2 SERPINE1",
              "M3_TwoCore": "M3 TwoCore", "M4_ThreeTop": "M4 ThreeTop",
              "M5_Panel10": "M5 Panel10"}


def f3(x):
    return "%.3f" % x


def ci(lo, hi, nd=3):
    return ("%.*f–%.*f" % (nd, lo, nd, hi)) if np.isfinite(lo) and np.isfinite(hi) else "—"


def pfmt(p):
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return "—"
    if p <= 0:
        return "<0.001"
    if p < 0.001:
        return "<0.001"
    return "%.3f" % p


rows = []            # (kind, c1, c2, c3, c4)
rows.append(("sec", "诊断模型判别性能（AUC；P 为 Hosmer–Lemeshow 校准检验）", "", "", ""))
for _, r in auc.iterrows():
    rows.append(("row", MODEL_NAME.get(r.model, r.model), f3(r.AUC),
                 ci(r.CI_low, r.CI_high), pfmt(r.HL_p)))

rows.append(("sec", "单因素 logistic 回归（每 1 个 log2 表达单位）", "", "", ""))
row = orr.loc["SERPINE1"]
rows.append(("row", "SERPINE1 OR", "%.2f" % row.OR,
             ci(row.OR_low, row.OR_high, 2), "%.3f" % row.p))

rows.append(("sec", "共线性诊断（方差膨胀因子 VIF）", "", "", ""))
for _, r in vif.iterrows():
    rows.append(("row", r.gene, "%.2f" % r.VIF, "—", "—"))

# ------------------------------------------------------------------ TSV
hdr = ["模型/指标", "值", "95%CI", "P"]
lines = ["\t".join(hdr)]
for kind, a, b, c, d in rows:
    lines.append("\t".join([a, b, c, d]))
tsv = os.path.join(OUT, "Table1_diagnostic.tsv")
with open(tsv, "w", encoding="utf-8-sig", newline="") as fh:
    fh.write("\n".join(lines) + "\n")
print("wrote", tsv)

# ------------------------------------------------------------------ preview
H_MT, H_HDR, H_ROW, H_SEC, H_FT, H_MB = 0.10, 0.30, 0.205, 0.255, 0.46, 0.06
H = H_MT + H_HDR + sum(H_SEC if r[0] == "sec" else H_ROW for r in rows) + H_FT + H_MB

fig = plt.figure(figsize=(6.3, H))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_axis_off()
ax.set_xlim(0, 1)
ax.set_ylim(0, H)

C1, C2, C3, C4 = 0.004, 0.545, 0.745, 0.996
y = H - H_MT
ax.plot([0, 1], [y, y], color="#111111", lw=1.15, solid_capstyle="butt")

# 表头
yh = y - H_HDR
ax.text(C1, yh + H_HDR * 0.45, hdr[0], fontsize=7.2, fontweight="bold",
        va="center", ha="left")
ax.text(C2, yh + H_HDR * 0.45, hdr[1], fontsize=7.2, fontweight="bold",
        va="center", ha="right")
ax.text(C3, yh + H_HDR * 0.45, hdr[2], fontsize=7.2, fontweight="bold",
        va="center", ha="center")
ax.text(C4, yh + H_HDR * 0.45, hdr[3], fontsize=7.2, fontweight="bold",
        va="center", ha="right")
ax.plot([0, 1], [yh, yh], color="#111111", lw=0.55, solid_capstyle="butt")
y = yh

for kind, a, b, c, d in rows:
    h = H_SEC if kind == "sec" else H_ROW
    ym = y - h / 2
    if kind == "sec":
        ax.add_patch(Rectangle((0, y - h), 1, h, facecolor="#F2F3F5",
                               edgecolor="none", zorder=0))
        ax.text(C1, ym, a, fontsize=6.9, fontweight="bold", va="center",
                ha="left", color="#1A1A1A")
    else:
        ax.text(C1, ym, a, fontsize=7.0, va="center", ha="left")
        ax.text(C2, ym, b, fontsize=7.0, va="center", ha="right")
        ax.text(C3, ym, c, fontsize=7.0, va="center", ha="center")
        ax.text(C4, ym, d, fontsize=7.0, va="center", ha="right")
    y -= h

ax.plot([0, 1], [y, y], color="#111111", lw=1.15, solid_capstyle="butt")

note = ("注：AUC 及其 95%CI 来自 ROC 分析，P 为 Hosmer–Lemeshow 校准检验；"
        "M1/M2 为单基因模型，M3–M5 分别为 2、3、10 基因组合。\n"
        "OR 为 SERPINE1 每升高 1 个 log2 表达单位的比值比；VIF 无对应显著性检验，"
        "VIF > 10 提示存在共线性。")
ax.text(0.004, y - 0.12, note, fontsize=5.7, va="top", ha="left",
        color="#444444", linespacing=1.55)

png = os.path.join(OUT, "Table1_preview.png")
fig.savefig(png, dpi=400, bbox_inches="tight", pad_inches=0.03,
            facecolor="white")
print("saved", png)
print("\n".join(lines))