# -*- coding: utf-8 -*-
"""图S3 五个基因虚拟敲除后的共享受扰靶点
a  共享敲除数分布（交集计数条）+ 10 个共享受扰基因 x 5 个敲除的 Z 矩阵
b  五个敲除的上皮程序扰动强度（18 上皮基因 mean |Z|）

注：含 $...$ 的字符串整体交给 mathtext，中文字形会缺字，
   故含数学符号的文本一律只用 ASCII，中文文本一律不含 $。
"""
import os

import numpy as np
import pandas as pd
import figstyle as S
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

tk = S.load(8, "D8_third_ko.tsv")
epi = S.load(8, "D8_epithelial_program_score.tsv")

TARGETS = ["Serpine1", "Tmsb4x", "P4hb", "S100a13", "Vim"]
SHOW = ["Krt5", "Perp", "Calml3", "S100a11", "Krt12", "Lypd2",
        "Gsta4", "Upk3bl", "Slurp1", "Anxa2"]

# 只保留非目标基因（排除敲除基因自身）
nt = tk[~tk["is_target"].astype(bool)].copy()
sig = nt[nt["p.value"] < 0.05]
sets = {t: set(sig[sig.target == t].gene) for t in TARGETS}
print("名义 P<0.05 的非目标显著数:", {t: len(sets[t]) for t in TARGETS})
print("FDR<0.05 的非目标显著数:",
      {t: int(((nt.target == t) & (nt["p.adj"] < 0.05)).sum()) for t in TARGETS})

deg = {g: sum(g in sets[t] for t in TARGETS) for g in nt.gene.unique()}
cnt = {d: sum(v == d for v in deg.values()) for d in range(1, 6)}
print("共享敲除数分布:", cnt)

zmap = dict(zip(zip(nt.gene, nt.target), nt.Z))
pmap = dict(zip(zip(nt.gene, nt.target), nt["p.value"]))

fig, axes = S.grid(1, 2, w=6.55, h=3.62, wspace=0.60, left=0.075, right=0.975,
                   top=0.885, bottom=0.155, widths=[1.46, 1.0])

# ---------------------------------------------------------------- (a)
axes[0].set_axis_off()
p = axes[0].get_position()
x0, y0, w0, h0 = p.x0, p.y0, p.width, p.height
ax_top = fig.add_axes([x0 + 0.010, y0 + h0 * 0.715, w0 * 0.80, h0 * 0.245])
ax_mat = fig.add_axes([x0 + 0.010, y0 + h0 * 0.035, w0 * 0.80, h0 * 0.545])

# 交集计数条：按共享敲除数分组
dvals = np.arange(1, 6)
bars = [cnt[d] for d in dvals]
ax_top.bar(dvals, bars, color=[S.PAL["lgrey"]] * 4 + [S.PAL["navy"]],
           width=0.62, edgecolor="none")
for d, v in zip(dvals, bars):
    ax_top.text(d, v + 0.9, str(v), ha="center", va="bottom", fontsize=5.8,
                color="#333333")
ax_top.set_xlim(0.35, 5.65)
ax_top.set_ylim(0, 29)
ax_top.set_xticks(dvals)
ax_top.set_xticklabels([f"{d}" for d in dvals], fontsize=6.0)
ax_top.set_yticks([0, 10, 20])
ax_top.set_xlabel("共享该靶点的敲除数", fontsize=6.2, labelpad=1.2)
ax_top.set_ylabel("基因数", fontsize=6.4, labelpad=2.5)
ax_top.tick_params(axis="both", labelsize=6.0, length=2.2)
S.ygrid(ax_top)
S.clean(ax_top)


# 基因 x 敲除 的 Z 矩阵
rows = SHOW
nR, nC = len(rows), len(TARGETS)
cmap = S.CMAP_SEQ
zmin, zmax = 0.9, 2.95
for i, g in enumerate(rows):
    yy = nR - 1 - i
    for j, t in enumerate(TARGETS):
        z = zmap.get((g, t), np.nan)
        pv = pmap.get((g, t), np.nan)
        on = (not np.isnan(pv)) and pv < 0.05
        if on:
            c = cmap((z - zmin) / (zmax - zmin))
        else:
            c = "#EFEFEF"
        ax_mat.add_patch(Rectangle((j - 0.46, yy - 0.42), 0.92, 0.84,
                                   facecolor=c, edgecolor="white", lw=0.8))
        txt = f"{z:.1f}"
        ax_mat.text(j, yy, txt, ha="center", va="center", fontsize=5.5,
                    color="white" if (on and z > 2.35) else "#2A2A2A")
    n = deg.get(g, 0)
    ax_mat.text(nC + 0.12, yy, f"{n}/5", ha="left", va="center", fontsize=5.8,
                color=S.PAL["navy"] if n >= 3 else "#666666",
                fontweight="bold" if n >= 3 else "normal")

ax_mat.set_xlim(-0.55, nC + 0.95)
ax_mat.set_ylim(-0.55, nR - 0.45)
ax_mat.set_xticks(range(nC))
ax_mat.set_xticklabels(TARGETS, fontsize=6.0)
ax_mat.set_yticks(range(nR))
ax_mat.set_yticklabels(rows[::-1], fontsize=6.0)
ax_mat.tick_params(length=0, pad=2)
for s in ax_mat.spines.values():
    s.set_visible(False)
ax_mat.text(-0.55, -1.62, "格内 = 虚拟敲除后该基因的 Z 值，右侧为该基因被几个敲除共享；"
            "灰底 = 名义 P ≥ 0.05；非目标 FDR < 0.05 者五个敲除均为 0",
            ha="left", va="top", fontsize=4.9, color="#666666")

# ---------------------------------------------------------------- (b)
ax = axes[1]
epi = epi.set_index("target").loc[TARGETS].reset_index()
y = np.arange(len(TARGETS))[::-1]
cols = [S.PAL["navy"], S.PAL["blue"], S.PAL["teal"], S.PAL["amber"],
        S.PAL["red"]]
ax.barh(y, epi["mean_Z_epi"], color=cols, height=0.62, edgecolor="none")
for yi, v, ns, ntst in zip(y, epi["mean_Z_epi"], epi["n_epi_sig"],
                           epi["n_epi_tested"]):
    ax.text(v + 0.045, yi, f"{v:.3f}", va="center", ha="left", fontsize=6.0)
    ax.text(2.46, yi, f"{int(ns)}/{int(ntst)}", va="center", ha="right",
            fontsize=5.4, color="#777777")
ax.text(2.46, 4.52, "显著/n", va="bottom", ha="right",
        fontsize=5.2, color="#999999")
ax.set_yticks(y)
ax.set_yticklabels(epi["target"], fontsize=6.8)
ax.set_xlim(0, 2.50)
ax.set_ylim(-0.55, 4.78)
ax.set_xticks([0, 0.5, 1.0, 1.5, 2.0])
ax.set_xlabel("18 个上皮基因的均值 |Z|", fontsize=6.6, labelpad=3)
ax.set_title("上皮程序扰动强度", pad=5, fontsize=7.6)
S.xgrid(ax) if hasattr(S, "xgrid") else None
ax.set_axisbelow(True)
ax.xaxis.grid(True, color="#EDEDED", lw=0.5)
S.clean(ax)
ax.text(-0.30, -0.135, "非目标显著 = 0（FDR < 0.05；仅靶基因自身达显著）",
        transform=ax.transAxes, ha="left", va="top", fontsize=5.2,
        color="#555555")
fig.text(0.0115, 0.925, "a", fontsize=11, fontweight="bold", va="bottom", ha="left", color="#111111")
fig.text(0.612, 0.925, "b", fontsize=11, fontweight="bold", va="bottom", ha="left", color="#111111")

S.save(fig, "S3_five_ko", also_pdf=False, tight=False)
# PDF 与 PNG 使用同一画布（避免 tight 裁切后上下留白与 PNG 不一致）
fig.savefig(os.path.join(S.OUT, "S3_five_ko.pdf"))
plt.close(fig)

