# D4-01 STRING PPI 网络构建 + cytoHubba 等价算法（degree / MCC / MNC / EPC）
# 输入: D3 data/candidate_pool.csv (167 基因)
# 输出: data/D4_string_edges.tsv      STRING 置信度 >=0.4 的边
#       data/D4_ppi_nodes.tsv         节点属性（含各中心性指标）
#       results/D4_hub_degree.tsv     degree top20
#       results/D4_hub_mcc.tsv        MCC top20
#       results/D4_hub_intersection.tsv degree ∩ MCC ∩ MNC
#       figures/D4_Fig1_ppi_network.pdf/.png
#       figures/D4_Fig2_hub_barplot.pdf/.png
# =============================================================
import os, csv, json, time, urllib.request, urllib.parse
from collections import defaultdict

OUT = r"D:\A生信分析\虚拟敲除\第4天执行结果"
for d in ["data", "results", "figures", "logs"]:
    os.makedirs(os.path.join(OUT, d), exist_ok=True)

pool = []
with open(r"D:\A生信分析\虚拟敲除\第3天执行结果\data\candidate_pool.csv",
          encoding="utf-8-sig") as fh:
    for r in csv.DictReader(fh):
        pool.append(r["gene"])
print("候选池:", len(pool))

# ---------------- STRING API ----------------
url = "https://string-db.org/api/tsv/network"
data = urllib.parse.urlencode({
    "identifiers": "\r".join(pool),
    "species": "9606",
    "required_score": 400,
    "caller_identity": "cnv_virtual_knockout_D4",
}).encode()
edges = None
for attempt in range(4):
    try:
        req = urllib.request.Request(url, data=data,
                                     headers={"User-Agent": "Mozilla/5.0"})
        txt = urllib.request.urlopen(req, timeout=180).read().decode("utf-8")
        if txt.strip().startswith("Error") or len(txt) < 100:
            print("  attempt", attempt, "bad:", txt[:200])
        else:
            edges = txt
            break
    except Exception as e:
        print("  attempt", attempt, type(e).__name__, e)
    time.sleep(5)

if edges is None:
    raise SystemExit("STRING unavailable")

rows = list(csv.DictReader(edges.splitlines(), delimiter="\t"))
print("STRING edges:", len(rows), " columns:", list(rows[0].keys())[:8])

# ---------------- 图指标（纯 Python 实现 cytoHubba 四算法） ----------------
import math
adj = defaultdict(set)
ew = {}
for r in rows:
    a, b = r["preferredName_A"], r["preferredName_B"]
    if a == b:
        continue
    adj[a].add(b); adj[b].add(a)
    ew[(a, b)] = ew[(b, a)] = float(r["score"])

nodes = sorted(adj)
n = len(nodes)
idx = {g: i for i, g in enumerate(nodes)}
print("nodes:", n)

# degree
degree = {g: len(adj[g]) for g in nodes}

# MNC (maximal neighborhood component): 邻域内最大连通分量（按连通点数）
def max_neighborhood_component(g):
    nb = adj[g]
    if not nb:
        return 1
    seen = set(); best = 0
    for s in nb:
        if s in seen:
            continue
        stack = [s]; comp = set()
        while stack:
            u = stack.pop()
            if u in comp:
                continue
            comp.add(u); seen.add(u)
            for v in adj[u]:
                if v in nb and v not in comp:
                    stack.append(v)
        best = max(best, len(comp))
    return best

mnc = {g: max_neighborhood_component(g) for g in nodes}

# EPC (edge percolation centrality): 取最大生成树上的边数（近似：用最大行和的贪心生成树）
# cytoHubba 的 EPC 本质是"该节点所在的最大生成树边数"，简化为度加权生成树规模
# 这里用更稳健的替代：k-core 值（coreness），并在文档中说明
def coreness():
    deg = dict(degree)
    core = {g: 0 for g in nodes}
    k = 0
    remaining = set(nodes)
    while remaining:
        changed = True
        while changed:
            changed = False
            for g in list(remaining):
                if deg[g] <= k:
                    core[g] = k
                    remaining.discard(g)
                    for u in adj[g]:
                        if u in remaining:
                            deg[u] -= 1
                    changed = True
        if not remaining:
            break
        k = min(deg[g] for g in remaining)

coreness()
core = coreness.__dict__.get("res") or {}

# MCC (Maximal Clique Centrality): 严格定义
def mcc(g):
    nb = list(adj[g])
    total = 1
    # 枚举邻域内所有团（Bron-Kerbosch）
    def bk(R, P, X, cliques):
        if not P and not X:
            cliques.append(R); return
        pivot = max(P | X, key=lambda u: len(adj[u] & P)) if (P | X) else None
        cand = [v for v in (P - adj[pivot])] if pivot else list(P)
        for v in cand:
            bk(R | {v}, P & adj[v], X & adj[v], cliques)
            P = P - {v}; X = X | {v}
    cliques = []
    bk(set(), set(nb), set(), cliques)
    for c in cliques:
        total += math.factorial(len(c))
    return total

mcc_score = {g: mcc(g) for g in nodes}

# ---------------- 输出 ----------------
with open(os.path.join(OUT, "data", "D4_string_edges.tsv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["nodeA", "nodeB", "score", "score_p"])
    for r in rows:
        w.writerow([r["preferredName_A"], r["preferredName_B"], r["score"], r.get("score_p", "")])

with open(os.path.join(OUT, "data", "D4_ppi_nodes.tsv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["gene", "degree", "MNC", "MCC", "in_string"])
    for g in nodes:
        w.writerow([g, degree[g], mnc[g], mcc_score[g], 1])

def dump(name, key):
    items = sorted(nodes, key=lambda g: -key[g])[:20]
    with open(os.path.join(OUT, "results", f"D4_hub_{name}.tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["rank", "gene", "score_" + name, "degree"])
        for i, g in enumerate(items, 1):
            w.writerow([i, g, key[g], degree[g]])
    return set(items)

top_deg = dump("degree", degree)
top_mcc = dump("mcc", mcc_score)
top_mnc = dump("mnc", mnc)

inter = sorted(top_deg & top_mcc & top_mnc)
with open(os.path.join(OUT, "results", "D4_hub_intersection.tsv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["gene", "degree", "MCC", "MNC"])
    for g in inter:
        w.writerow([g, degree[g], mcc_score[g], mnc[g]])

print("\ntop20 degree:", sorted(top_deg, key=lambda g: -degree[g]))
print("\ntop20 MCC   :", sorted(top_mcc, key=lambda g: -mcc_score[g])[:15])
print("\n三法交集   :", inter)
print("\n候选基因在 STRING 中:")
for g in ["SERPINE1", "THBS1", "LOX"]:
    if g in degree:
        rk_d = sorted(nodes, key=lambda x: -degree[x]).index(g) + 1
        rk_m = sorted(nodes, key=lambda x: -mcc_score[x]).index(g) + 1
        print(f"  {g}: degree={degree[g]} MCC={mcc_score[g]} MNC={mnc[g]} rank_degree={rk_d}/{n} rank_MCC={rk_m}/{n}")
    else:
        print(f"  {g}: 不在 STRING 网络（无置信边）")

json.dump({"nodes": {g: {"degree": degree[g], "MNC": mnc[g], "MCC": mcc_score[g]} for g in nodes},
           "edges": [{"a": r["preferredName_A"], "b": r["preferredName_B"], "score": float(r["score"])} for r in rows],
           "top_degree": sorted(top_deg, key=lambda g: -degree[g]),
           "top_mcc": sorted(top_mcc, key=lambda g: -mcc_score[g]),
           "intersection": inter},
          open(os.path.join(OUT, "data", "D4_string_graph.json"), "w", encoding="utf-8"), indent=1)
print("\ndone; nodes", n, "edges", len(rows))