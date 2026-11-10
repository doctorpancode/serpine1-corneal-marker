# -*- coding: utf-8 -*-
"""D4-02 PPI 网络图 + Hub 条形图（matplotlib，无 igraph/networkx 依赖的力导向布局）
输入: data/D4_string_edges.tsv, data/D4_ppi_nodes.tsv, data/D4_string_graph.json
输出: figures/D4_Fig1_ppi_network.pdf/.png, figures/D4_Fig2_hub_barplot.pdf/.png
"""
import os, json, csv, math, random
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle
import numpy as np

OUT = r"D:\A生信分析\虚拟敲除\第4天执行结果"
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42

G = json.load(open(os.path.join(OUT, "data", "D4_string_graph.json"), encoding="utf-8"))
nodes = list(G["nodes"].keys())
edges = G["edges"]
ni = {g: i for i, g in enumerate(nodes)}
deg = np.array([G["nodes"][g]["degree"] for g in nodes])
mcc = np.array([G["nodes"][g]["MCC"] for g in nodes])

# ---- 力导向布局（Fruchterman-Reingold，固定种子） ----
random.seed(42); np.random.seed(42)
n = len(nodes)
pos = np.random.rand(n, 2) * 10
adj = [[] for _ in range(n)]
for e in edges:
    if e["a"] in ni and e["b"] in ni:
        a, b = ni[e["a"]], ni[e["b"]]
        adj[a].append(b); adj[b].append(a)

k = 1.4
for it in range(320):
    temp = 1.2 * (1 - it / 320) + 0.02
    disp = np.zeros((n, 2))
    diff = pos[:, None, :] - pos[None, :, :]
    dist = np.sqrt((diff ** 2).sum(-1)) + 1e-6
    np.fill_diagonal(dist, np.inf)
    rep = (k ** 2 / dist ** 2)[:, :, None] * diff
    disp += rep.sum(1)
    for a in range(n):
        for b in adj[a]:
            d = pos[a] - pos[b]
            dd = np.linalg.norm(d) + 1e-6
            disp[a] -= (dd ** 2 / k) * d / dd
    ln = np.linalg.norm(disp, axis=1) + 1e-9
    pos += disp / ln[:, None] * np.minimum(ln, temp)[:, None]
pos -= pos.mean(0)
pos /= np.abs(pos).max()

inter = set(G["intersection"])
sub = [g for g in ["SERPINE1", "THBS1", "LOX"] if g in ni]

fig, ax = plt.subplots(figsize=(11, 10))
segs, ws = [], []
for e in edges:
    if e["a"] in ni and e["b"] in ni:
        segs.append([pos[ni[e["a"]]], pos[ni[e["b"]]]])
        ws.append(e["score"])
lc = LineCollection(segs, array=np.array(ws), cmap="Greys", linewidths=0.35 + 1.6 * (np.array(ws) - 0.4) / 0.6, alpha=0.55)
ax.add_collection(lc)
cb = fig.colorbar(lc, ax=ax, fraction=0.025, pad=0.01)
cb.set_label("STRING combined score", fontsize=9)

sz = 60 + 260 * (deg - deg.min()) / (deg.max() - deg.min())
cols = ["#d62728" if g in inter else ("#ff7f0e" if g in sub else "#8fb8de") for g in nodes]
ax.scatter(pos[:, 0], pos[:, 1], s=sz, c=cols, edgecolors="white", linewidths=0.8, zorder=3, alpha=0.92)

lab = [g for g in nodes if g in inter or g in sub or deg[ni[g]] >= np.percentile(deg, 97)]
for g in lab:
    ax.annotate(g, pos[ni[g]], fontsize=7.5, fontweight="bold", ha="center", va="center",
                zorder=4, xytext=(0, 9), textcoords="offset points")
ax.set_title("STRING PPI network of the D3 candidate pool (167 genes -> %d nodes / %d edges)\n"
             "red = degree AND MCC AND MNC intersection; orange = prior candidates"
             % (n, len(edges)), fontsize=11)
ax.set_xlim(pos[:, 0].min() - 0.12, pos[:, 0].max() + 0.12)
ax.set_ylim(pos[:, 1].min() - 0.12, pos[:, 1].max() + 0.12)
ax.axis("off")
for ext in ["pdf", "png"]:
    fig.savefig(os.path.join(OUT, "figures", "D4_Fig1_ppi_network." + ext),
                dpi=300, bbox_inches="tight")
plt.close(fig)

# ---- Fig2: degree / MCC 双条形图 ----
order = sorted(nodes, key=lambda g: -deg[ni[g]])[:20]
fig, axes = plt.subplots(1, 2, figsize=(13, 6.2), gridspec_kw={"width_ratios": [1, 1]})
for ax, key, val, ttl in [
    (axes[0], "Node degree", [G["nodes"][g]["degree"] for g in order], "cytoHubba-style degree (top 20)"),
    (axes[1], "MCC score (log10)", [math.log10(G["nodes"][g]["MCC"]) for g in order], "Maximal Clique Centrality (same genes)"),
]:
    c = ["#d62728" if g in inter else "#4c72b0" for g in order]
    ax.barh(range(len(order))[::-1], val, color=c, edgecolor="white")
    ax.set_yticks(range(len(order))[::-1])
    ax.set_yticklabels(order, fontsize=9)
    ax.set_xlabel(key, fontsize=10)
    ax.set_title(ttl, fontsize=11)
    ax.grid(axis="x", alpha=0.25)
fig.suptitle("Hub gene ranking of the D3 candidate pool (red = three-method intersection)", fontsize=12)
fig.tight_layout()
for ext in ["pdf", "png"]:
    fig.savefig(os.path.join(OUT, "figures", "D4_Fig2_hub_barplot." + ext), dpi=300, bbox_inches="tight")
plt.close(fig)
print("figures done; intersec =", sorted(inter))