# -*- coding: utf-8 -*-
"""Figure S9: transport calibration + decision curve for SERPINE1 (GSE58291, n = 27)."""
import io
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
R = r"D:\A生信分析\虚拟敲除"
OUT = R + r"\文章\figs_v2\out\FigS9_calibration.png"

plt.rcParams.update({
    "font.family": "Arial", "font.size": 8, "axes.linewidth": 0.7,
    "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "axes.spines.top": False, "axes.spines.right": False,
})
BLUE, ORANGE, GREY, RED = "#3B6FB6", "#E08214", "#8C8C8C", "#C0392B"


def logit_fit(X, y, iters=300):
    b = np.zeros(X.shape[1])
    for _ in range(iters):
        eta = np.clip(X @ b, -30, 30)
        p = 1 / (1 + np.exp(-eta))
        W = np.clip(p * (1 - p), 1e-9, None)
        H = X.T @ (X * W[:, None])
        gr = X.T @ (y - p)
        try:
            step = np.linalg.solve(H, gr)
        except np.linalg.LinAlgError:
            break
        b = b + step
        if np.max(np.abs(step)) < 1e-9:
            break
    return b


# discovery model on z-scored SERPINE1
expr = pd.read_csv(R + r"\D1执行结果\第一天执行结果\data\processed\GSE241715_expression.tsv", sep="\t").set_index("SYMBOL")
sheet = pd.read_csv(R + r"\D1执行结果\第一天执行结果\data\processed\GSE241715_sample_sheet.tsv", sep="\t")
y1 = sheet.set_index("sample").loc[expr.columns, "binary_lesion"].values.astype(int)
x1 = expr.loc["SERPINE1"].values.astype(float)
z1 = (x1 - x1.mean()) / x1.std(ddof=1)
b1 = logit_fit(np.column_stack([np.ones_like(z1), z1]), y1)

g58 = pd.read_csv(R + r"\第9天执行结果\data\GSE58291_expr.tsv", sep="\t", low_memory=False).set_index("SYMBOL")
grp = pd.read_csv(R + r"\第9天执行结果\results\D9_GSE58291_sample_group.tsv", sep="\t")
g58 = g58[[c for c in g58.columns if str(c).endswith("AVG_Signal")]]
g58.columns = [re.sub(r"\s+", "", str(c).split(".")[0]) for c in g58.columns]
gmap = dict(zip(grp.iloc[:, 0].astype(str), grp.iloc[:, 1].astype(str)))
keep = [c for c in g58.columns if c in gmap]
x2 = g58.loc["SERPINE1", keep].values.astype(float)
y2 = np.array([0 if gmap[c].lower().startswith("cont") else 1 for c in keep])
z2 = (x2 - x2.mean()) / x2.std(ddof=1)
p2 = 1 / (1 + np.exp(-(b1[0] + b1[1] * z2)))

fig = plt.figure(figsize=(6.6, 2.35))
gs = GridSpec(1, 3, width_ratios=[1, 1, 1], wspace=0.38, left=0.075, right=0.985,
              bottom=0.205, top=0.9)

# (a) calibration by tertiles of predicted risk
ax = fig.add_subplot(gs[0])
q = pd.qcut(p2, 3, labels=False, duplicates="drop")
xs, ys, los, his, ns = [], [], [], [], []
for g in np.unique(q):
    m = q == g
    xs.append(p2[m].mean())
    k, n = int(y2[m].sum()), int(m.sum())
    ys.append(k / n)
    lo, hi = 0.0, 1.0
    # Wilson interval
    zc = 1.96
    den = 1 + zc ** 2 / n
    centre = (k / n + zc ** 2 / (2 * n)) / den
    half = zc * np.sqrt(k / n * (1 - k / n) / n + zc ** 2 / (4 * n ** 2)) / den
    los.append(max(0.0, centre - half))
    his.append(min(1.0, centre + half))
    ns.append(n)
ax.plot([0, 1], [0, 1], ls="--", lw=0.8, color=GREY, label="Ideal")
ax.errorbar(xs, ys, yerr=[np.array(ys) - los, np.array(his) - np.array(ys)], fmt="o-",
            color=BLUE, ms=4, lw=1.2, capsize=2.5, label="Observed (tertiles)")
slope, inter = 1.191, -0.757
grid = np.linspace(0.02, 0.98, 50)
ax.plot(grid, 1 / (1 + np.exp(-(inter + slope * np.log(grid / (1 - grid))))), color=ORANGE, lw=1,
        label="Slope 1.19, intercept −0.76")
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_xlabel("Predicted probability")
ax.set_ylabel("Observed proportion")
ax.set_title("a  Calibration (n = 27)", loc="left", fontsize=8.5, fontweight="bold")
ax.legend(frameon=False, fontsize=6.2, loc="upper left")
ax.text(0.03, 0.62, "HL χ² = 10.79\nP = 0.010", fontsize=6.4, color=GREY)

# (b) decision curve
ax = fig.add_subplot(gs[1])
prev = y2.mean()
th = np.arange(0.05, 0.96, 0.01)
nb, nb_all = [], []
for t in th:
    flag = p2 >= t
    tp = np.sum(flag & (y2 == 1))
    fp = np.sum(flag & (y2 == 0))
    nb.append(tp / len(y2) - fp / len(y2) * t / (1 - t))
    nb_all.append(prev - (1 - prev) * t / (1 - t))
ax.plot(th, nb, color=BLUE, lw=1.3, label="SERPINE1 model")
ax.plot(th, nb_all, color=GREY, ls="--", lw=0.9, label="Treat all")
ax.axhline(0, color=RED, ls=":", lw=0.9, label="Treat none")
ax.set_ylim(-0.35, 0.6)
ax.set_xlim(0.05, 0.95)
ax.set_xlabel("Threshold probability")
ax.set_ylabel("Net benefit")
ax.set_title("b  Decision curve", loc="left", fontsize=8.5, fontweight="bold")
ax.legend(frameon=False, fontsize=6.2, loc="lower left")
ax.annotate("NB 0.506\n@ 0.40", xy=(0.40, 0.506), xytext=(0.5, 0.52), fontsize=6.2, color=BLUE,
            arrowprops=dict(arrowstyle="-", lw=0.6, color=BLUE))

# (c) predicted risk by group
ax = fig.add_subplot(gs[2])
bin_ = np.arange(0, 1.01, 0.1)
for i, (lab, col, name) in enumerate([(0, GREY, "Control"), (1, ORANGE, "Infected")]):
    v = p2[y2 == lab]
    ax.hist(v, bins=bin_, color=col, alpha=0.75, label=f"{name} (n = {len(v)})",
            edgecolor="white", linewidth=0.6)
ax.set_xlim(0, 1)
ax.set_xlabel("Predicted probability")
ax.set_ylabel("Samples")
ax.set_title("c  Risk distribution", loc="left", fontsize=8.5, fontweight="bold")
ax.legend(frameon=False, fontsize=6.2, loc="upper center")

fig.savefig(OUT, dpi=400)
fig.savefig(OUT.replace(".png", ".pdf"))
print("saved:", OUT)
print("tertiles:", [(round(float(a), 3), round(float(b), 3), c) for a, b, c in zip(xs, ys, ns)])
print("mean pred", round(float(p2.mean()), 3), "obs", round(float(y2.mean()), 3))