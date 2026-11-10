# D3-01 候选池构建：主题分层三源并集 → 与 WGCNA 模块收口
# 输入:
#   work\ctd_cornea_sets.json   241 个疾病基因集 (CTD 34 + DisGeNET 89 + HPO 51 + DISEASES2025 68)
#   D2 data\wgcna_significant_module_genes.tsv
#   D1 GSE241715_expression.rds（用于与差异基因取交）
# 输出:
#   results\D3_set_inventory.tsv         所有用到的基因集规模
#   results\D3_theme_A_CNV.tsv          主题A：CNV 直接证据
#   results\D3_theme_B_pathological.tsv 主题B：病理性新生血管泛化 × 角膜语境收口
#   results\D3_theme_C_corneal_injury.tsv 主题C：角膜/角膜缘损伤 × 模块收口
#   results\D3_venn_counts.tsv          三主题两两/三交规模
#   data\candidate_pool.csv             最终候选池（含来源标记）
#   figures\D3_Fig1_venn.pdf/.png
# =============================================================
import json, os, csv, re
from collections import Counter

WORK = r"D:\A生信分析\虚拟敲除\work"
D2   = r"D:\A生信分析\虚拟敲除\第2天执行结果"
D1   = r"D:\A生信分析\虚拟敲除\D1执行结果\D1_CNV_角膜新生血管_数据分析"
OUT  = r"D:\A生信分析\虚拟敲除\第3天执行结果"
for d in ["data", "results", "figures", "logs"]:
    os.makedirs(os.path.join(OUT, d), exist_ok=True)

sets = json.load(open(os.path.join(WORK, "ctd_cornea_sets.json"), encoding="utf-8"))
print("总基因集:", len(sets))

def U(keys):
    s = set()
    for k in keys:
        s |= set(sets[k])
    return s

# ---------- 主题 A：CNV 直接证据 ----------
A_keys = [k for k in sets if re.search(r"neovascular", k, re.I)
          and not re.search(r"choroidal|retinal|vitreoretinopathy|glaucoma|myopic|optic disc", k, re.I)]
print("\n主题A 使用集:")
for k in sorted(A_keys):
    print(f"  {len(sets[k]):5d}  {k}")
A = U(A_keys)
print("主题A 并集:", len(A))

# ---------- 主题 B：病理性新生血管泛化 × 角膜语境 ----------
B_raw_keys = [k for k in sets if re.search(r"^Neovascularization, Pathologic", k)]
B_raw = U(B_raw_keys)
# 角膜语境收口：与"所有角膜相关集"的并集取交，去掉视网膜/肿瘤语境基因
cornea_keys = [k for k in sets if re.search(r"cornea|corneal|keratitis|keratoplasty|limbal", k, re.I)]
cornea_union = U(cornea_keys)
B = B_raw & cornea_union
print(f"\n主题B: Neovascularization,Pathologic = {len(B_raw)} ; 角膜语境内 = {len(B)}")
print("  角膜相关集数:", len(cornea_keys))

# ---------- 主题 C：角膜/角膜缘损伤 × WGCNA 模块 ----------
C_raw_keys = [k for k in sets if re.search(r"cornea|corneal|keratitis|keratoplasty|limbal|corneal ulcer|corneal disease",
                                           k, re.I)
              and not re.search(r"dystroph|Schnyder|Meesmann|Fleck|Thiel|Lattice|Granular|Avellino|gelatinous|macular|Sclerocornea|plana|ectasia|staphyloma", k, re.I)]
C_raw = U(C_raw_keys)
mod = {}
with open(os.path.join(D2, "data", "wgcna_significant_module_genes.tsv"), encoding="utf-8") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        mod.setdefault(r["gene"], set()).add(r["module"])
mod_genes = set(mod)
C = C_raw & mod_genes
print(f"\n主题C: 角膜/角膜缘损伤集 = {len(C_raw)} ; ∩ 显著模块基因({len(mod_genes)}) = {len(C)}")

# ---------- 三主题关系 ----------
venn = {
    "A": len(A), "B": len(B), "C": len(C),
    "A&B": len(A & B), "A&C": len(A & C), "B&C": len(B & C),
    "A&B&C": len(A & B & C),
    "A|B|C": len(A | B | C),
}
print("\nVenn:", venn)

# ---------- 候选池 ----------
pool = A | B | C
print("候选池（并集）:", len(pool))

src = {}
for g in pool:
    tags = []
    if g in A: tags.append("A:CNV_direct")
    if g in B: tags.append("B:pathological_neo")
    if g in C: tags.append("C:corneal_injury_x_module")
    src[g] = tags

# 每个基因的支持集数量（来自几个原始基因集）
sup = Counter()
for k, v in sets.items():
    for g in set(v) & pool:
        sup[g] += 1

# 表达式过滤：是否在 GSE241715 检出
genes_expr = set()
with open(os.path.join(D1, "data", "processed", "GSE241715_expression.tsv"), encoding="utf-8") as fh:
    hdr = fh.readline()
    for line in fh:
        genes_expr.add(line.split("\t")[0].strip('"'))
print("GSE241715 基因数:", len(genes_expr))

rows = []
for g in sorted(pool):
    rows.append({
        "gene": g,
        "n_source_sets": sup[g],
        "themes": "|".join(src[g]),
        "in_theme_A": int(g in A),
        "in_theme_B": int(g in B),
        "in_theme_C": int(g in C),
        "wgcna_modules": "|".join(sorted(mod.get(g, []))),
        "in_GSE241715": int(g in genes_expr),
    })
rows.sort(key=lambda r: (-r["n_source_sets"], -len(r["themes"]), r["gene"]))

with open(os.path.join(OUT, "data", "candidate_pool.csv"), "w", newline="", encoding="utf-8-sig") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)

# ---------- 明细表 ----------
with open(os.path.join(OUT, "results", "D3_set_inventory.tsv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["theme", "set_name", "n_genes", "n_after_theme"])
    for k in sorted(A_keys): w.writerow(["A", k, len(sets[k]), ""])
    for k in sorted(B_raw_keys): w.writerow(["B", k, len(sets[k]), len(B)])
    for k in sorted(C_raw_keys): w.writerow(["C", k, len(sets[k]), len(C)])

for nm, s in [("A", A), ("B", B), ("C", C)]:
    with open(os.path.join(OUT, "results", f"D3_theme_{nm}.tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene"])
        for g in sorted(s): w.writerow([g])

with open(os.path.join(OUT, "results", "D3_venn_counts.tsv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["region", "n_genes"])
    for k, v in venn.items(): w.writerow([k, v])

json.dump({"A": sorted(A), "B": sorted(B), "C": sorted(C), "pool": sorted(pool)},
          open(os.path.join(OUT, "data", "D3_sets.json"), "w", encoding="utf-8"), indent=1)

# 候选基因是否在池中 + 支持度
for g in ["SERPINE1", "THBS1", "LOX"]:
    print(f"  {g}: in_pool={g in pool} n_sets={sup.get(g,0)} themes={src.get(g)}")

print("\n完成。候选池:", len(pool), " 文件:", os.path.join(OUT, "data", "candidate_pool.csv"))