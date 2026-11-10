# -*- coding: utf-8 -*-
"""统一绘图样式与数据接口（角膜新生血管多组学稿件配图）。
所有图件由本模块统一控制版式：白底、无顶/右边框、7pt 刻度、8pt 轴标题、
8.5pt 面板标题、粗体面板字母；PNG 400dpi + PDF 矢量双导出。
"""
from __future__ import annotations

import os
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
import matplotlib.patheffects as pe

# ---------------------------------------------------------------- paths
_HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.environ.get("VKO_BASE") or os.path.abspath(os.path.join(_HERE, "..", ".."))
OUT = os.environ.get("VKO_OUT") or os.path.join(BASE, "figures_out")
os.makedirs(OUT, exist_ok=True)

META = os.path.join(BASE, "data", "metadata")
RAW = os.path.join(BASE, "data", "raw")

D2 = os.path.join(BASE, "results", "D2")
D3 = os.path.join(BASE, "results", "D3")
D4 = os.path.join(BASE, "results", "D4")
D5 = os.path.join(BASE, "results", "D5")
D6 = os.path.join(BASE, "results", "D6")
D7 = os.path.join(BASE, "results", "D7")
D8 = os.path.join(BASE, "results", "D8")
D9 = os.path.join(BASE, "results", "D9")

FIGDIR = OUT


def asset(name):
    """Locate a derived input: metadata -> raw -> figures_out."""
    for root in (META, RAW, OUT):
        cand = os.path.join(root, name)
        if os.path.exists(cand):
            return cand
    raise FileNotFoundError(name)


def asset_raw(name):
    """Locate a downloaded GEO matrix under data/raw."""
    for root in (RAW, META):
        cand = os.path.join(root, name)
        if os.path.exists(cand):
            return cand
    raise FileNotFoundError(name)


def load(day: int, fname: str, **kw) -> pd.DataFrame:
    p = os.path.join({
        2: D2, 3: D3, 4: D4, 5: D5, 6: D6, 7: D7, 8: D8, 9: D9
    }[day], fname)
    return pd.read_csv(p, sep=kw.pop("sep", "\t"), **kw)


# ---------------------------------------------------------------- style
PAL = dict(
    blue="#3B6BA5", red="#C0504D", grey="#8C8C8C", lgrey="#D9D9D9",
    teal="#2A9D8F", amber="#E9A23B", purple="#7B5EA7", green="#5A8F5A",
    navy="#2F4B7C", brick="#B4483F", sand="#C9A227",
)
GROUP_COLORS = {
    "Control": PAL["blue"], "uninjured": PAL["blue"], "uninjured ": PAL["blue"],
    "Viral": "#7FA8D4", "Keratoconus": "#B9CBE0",
    "Bacterial_fungal": PAL["red"], "NaOH": PAL["red"],
}
# 四个 bulk 组按疾病严重度排序
BULK_GROUPS = ["Control", "Keratoconus", "Viral", "Bacterial_fungal"]
BULK_LABEL = {"Control": "正常", "Keratoconus": "圆锥角膜",
              "Viral": "病毒性", "Bacterial_fungal": "细菌/真菌性"}
BULK_COLORS = ["#4C72B0", "#8FB0D0", "#E0A16B", "#C0504D"]

CMAP_DIV = LinearSegmentedColormap.from_list(
    "div", ["#2F4B7C", "#7BA6D0", "#F7F7F7", "#E39A93", "#9E2B25"])
CMAP_SEQ = LinearSegmentedColormap.from_list(
    "seq", ["#F7F8FA", "#C9DCEC", "#7FA8C9", "#3B6BA5", "#20375A"])

rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 400,
    "font.family": ["Arial", "Microsoft YaHei"],
    "font.size": 7.5,
    "axes.titlesize": 8.5,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "axes.linewidth": 0.7,
    "axes.edgecolor": "#333333",
    "axes.labelcolor": "#1A1A1A",
    "text.color": "#1A1A1A",
    "xtick.color": "#333333",
    "ytick.color": "#333333",
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "xtick.major.size": 2.6,
    "ytick.major.size": 2.6,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "legend.frameon": False,
    "legend.handlelength": 1.1,
    "legend.handletextpad": 0.5,
    "legend.borderpad": 0.2,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
    "mathtext.default": "regular",
    "lines.linewidth": 1.0,
    "patch.linewidth": 0.6,
})


# ---------------------------------------------------------------- helpers
def fig_ax(w=6.3, h=4.0, **kw):
    fig = plt.figure(figsize=(w, h))
    ax = fig.add_axes(kw.pop("rect", [0.12, 0.15, 0.84, 0.78]), **kw)
    return fig, ax


def grid(nrows, ncols, w=6.3, h=4.0, wspace=0.45, hspace=0.5,
         left=0.09, right=0.975, top=0.94, bottom=0.10, heights=None, widths=None):
    fig, axes = plt.subplots(nrows, ncols, figsize=(w, h),
                             gridspec_kw=dict(wspace=wspace, hspace=hspace,
                                              height_ratios=heights,
                                              width_ratios=widths))
    fig.subplots_adjust(left=left, right=right, top=top, bottom=bottom)
    return fig, np.asarray(axes).ravel()


def panel(ax, letter, x=-0.16, y=1.06, size=11):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=size,
            fontweight="bold", va="bottom", ha="left", color="#111111")


def clean(ax, left=True, bottom=True):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(left)
    ax.spines["bottom"].set_visible(bottom)
    return ax


def ygrid(ax, color="#EDEDED", lw=0.5):
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color=color, lw=lw)
    ax.xaxis.grid(False)


def xgrid(ax, color="#EDEDED", lw=0.5):
    ax.set_axisbelow(True)
    ax.xaxis.grid(True, color=color, lw=lw)
    ax.yaxis.grid(False)


def stars(p):
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return ""
    if p < 1e-4:
        return "****"
    if p < 1e-3:
        return "***"
    if p < 1e-2:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def fnum(x, nd=2):
    """科学计数友好格式。"""
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "NA"
    if x == 0:
        return "0"
    a = abs(x)
    if a >= 1e4 or a < 1e-3:
        e = int(np.floor(np.log10(a)))
        m = x / 10 ** e
        return f"{m:.1f}×10$^{{{e}}}$"
    if a >= 100:
        return f"{x:.0f}"
    if a >= 1:
        return f"{x:.2f}"
    return f"{x:.3f}"


def save(fig, name, also_pdf=True, tight=True):
    png = os.path.join(OUT, f"{name}.png")
    if tight:
        fig.savefig(png, bbox_inches="tight", pad_inches=0.02)
    else:
        fig.savefig(png)
    if also_pdf:
        fig.savefig(os.path.join(OUT, f"{name}.pdf"), bbox_inches="tight",
                    pad_inches=0.02)
    plt.close(fig)
    print("saved", png)
    return png


def hbar(ax, labels, values, colors, err=None, xlabel="", fmt="{:.2f}", pad=0.02):
    y = np.arange(len(labels))[::-1]
    ax.barh(y, values, color=colors, height=0.68, edgecolor="none",
            xerr=err, error_kw=dict(lw=0.6, capsize=1.6, ecolor="#555555"))
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel(xlabel)
    vmax = np.nanmax(values)
    for yi, v in zip(y, values):
        ax.text(v + vmax * pad, yi, fmt.format(v), va="center", ha="left", fontsize=6.3)
    ax.set_xlim(0, vmax * 1.18)
    clean(ax)
    return ax


def heat(ax, mat, rowlabels, collabels, cmap=CMAP_DIV, vmin=None, vmax=None,
         center=0, cbar_label="", cbar_kw=None, annot=None, annot_size=6,
         xrot=45):
    if center is not None and vmin is None and vmax is None:
        lim = np.nanmax(np.abs(mat))
        norm = TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim)
    else:
        norm = None
    im = ax.imshow(mat, cmap=cmap, norm=norm, vmin=None if norm else vmin,
                   vmax=None if norm else vmax, aspect="auto")
    ax.set_xticks(range(len(collabels)))
    ax.set_xticklabels(collabels, rotation=xrot, ha="right")
    ax.set_yticks(range(len(rowlabels)))
    ax.set_yticklabels(rowlabels)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    if annot is not None:
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                if annot[i][j] not in ("", None):
                    ax.text(j, i, annot[i][j], ha="center", va="center",
                            fontsize=annot_size, color="#1A1A1A")
    kw = dict(fraction=0.045, pad=0.02)
    kw.update(cbar_kw or {})
    cb = plt.colorbar(im, ax=ax, **kw)
    cb.outline.set_linewidth(0.5)
    cb.ax.tick_params(labelsize=6.2, width=0.5, length=2)
    if cbar_label:
        cb.set_label(cbar_label, fontsize=6.8)
    return im


def slopegraph(ax, names, v1, v2, l1, l2, c1=None, c2=None, fmt="{:.0f}"):
    y = np.arange(len(names))[::-1]
    for yi, a, b in zip(y, v1, v2):
        col = PAL["red"] if b < a else PAL["blue"]
        ax.plot([0, 1], [a, b], color=col, lw=1.0, marker="o", ms=3.0, zorder=3)
        ax.text(-0.04, a, fmt.format(a), ha="right", va="center", fontsize=6)
        ax.text(1.04, b, fmt.format(b), ha="left", va="center", fontsize=6)
    ax.set_xlim(-0.35, 1.35)
    ax.set_ylim(min(min(v1), min(v2)) - 20, max(max(v1), max(v2)) + 20)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([l1, l2])
    ax.set_yticks(y)
    ax.set_yticklabels(names)
    clean(ax)
    ax.spines["bottom"].set_visible(False)