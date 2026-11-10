# ============================================================
# D9-06  待办 10.3-(5)：换真 CellChat（替代 D6 的自实现 L-R 打分）
#
#   D6-06 用手写 LR 表 + 均值乘积当通讯强度，缺 CellChat 的核心功能：
#   - 数据库（小鼠 CellChatDB，2029 个互作）
#   - 通讯概率模型（Hill function + 质量作用定律）
#   - 置换检验（permutation test）给出 p 值
#   - 通路级聚合（信号通路而不是单个 LR 对）
#
#   本脚本在 D6 的 Seurat 对象上跑正式 CellChat，比较 NaOH vs uninjury，
#   并输出 "炎症细胞 → 血管细胞" 的主要信号通路（这是 Fig9 机制图的依据）。
#
#   注意：CellChat 2.2.0.9001 装自 GitHub main（CRAN 归档版已下架）。
#         数据库对象在 2.2 里是 list，用 CellChatDB.mouse$interaction。
#
#   产出：
#     data/D9_cellchat_NaOH.rds, data/D9_cellchat_uninj.rds
#     results/D9_cellchat_compare_pathways.tsv
#     results/D9_cellchat_inflam_to_vessel_real.tsv
#     results/D9_cellchat_lr_to_vessel_real.tsv
#     figures/main/D9_Fig9_cellchat_*.png/pdf
# ============================================================
suppressPackageStartupMessages({
  library(Seurat); library(CellChat); library(data.table)
  library(ggplot2); library(patchwork); library(ComplexHeatmap)
})
set.seed(20261110)
OUT <- "D:/A生信分析/虚拟敲除/第9天执行结果"
RES <- file.path(OUT, "results"); FIG <- file.path(OUT, "figures", "main")
DAT <- file.path(OUT, "data")

sc <- readRDS("D:/A生信分析/虚拟敲除/第6天执行结果/data/D6_seurat_GSE275678_final.rds")
sc <- JoinLayers(sc)
cat("Seurat:", dim(sc), "\n")
cat("细胞类型:\n"); print(table(sc$celltype))

# CellChat 需要 data 层（normalized），D6 对象已有
data.input <- GetAssayData(sc, layer = "data")
meta <- data.frame(labels = as.character(sc$celltype),
                   samples = as.character(sc$sample),
                   row.names = colnames(sc))
cat("分组:\n"); print(table(meta$samples))
rm(sc); gc()

db <- CellChatDB.mouse
cat("数据库互作数:", nrow(db$interaction), "\n")

run_cc <- function(grp) {
  cells <- rownames(meta)[meta$samples == grp]
  cc <- createCellChat(object = data.input[, cells], meta = meta[cells, , drop = FALSE],
                       group.by = "labels")
  cc@DB <- db
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc, thresh.p = 0.05, do.fast = TRUE)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = "triMean", raw.use = TRUE, population.size = TRUE)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc
}

cc_naoh <- run_cc("NaOH")
cat("[done] NaOH: interactions =", nrow(cc_naoh@net$prob |> as.data.frame()), "\n")
saveRDS(cc_naoh, file.path(DAT, "D9_cellchat_NaOH.rds"))

cc_un <- run_cc("uninjury")
cat("[done] uninjury\n")
saveRDS(cc_un, file.path(DAT, "D9_cellchat_uninjury.rds"))

## ---------- 1. 通路级比较 ----------
# CellChat 2.2 的 netP$prob 是 [source, target, pathway]（11 x 11 x 94），
# 不是 [pathway, source, target]。之前的索引写反导致 strength 实际是"细胞类型的总量"。
pw <- function(cc) {
  P <- cc@netP$prob                       # source x target x pathway
  pws <- dimnames(P)[[3]]
  rbindlist(lapply(pws, function(k) data.table(pathway = k, strength = sum(P[, , k], na.rm = TRUE))))
}
net <- merge(pw(cc_naoh), pw(cc_un),
             by = "pathway", suffixes = c("_NaOH", "_uninjured"), all = TRUE)
net[is.na(strength_NaOH), strength_NaOH := 0]
net[is.na(strength_uninjured), strength_uninjured := 0]
net[, log2FC := log2((strength_NaOH + 1e-6) / (strength_uninjured + 1e-6))]
net[, ratio := (strength_NaOH + 1e-6) / (strength_uninjured + 1e-6)]
setorder(net, -strength_NaOH)
fwrite(net, file.path(RES, "D9_cellchat_compare_pathways.tsv"), sep = "\t")
cat("\n=== 通路强度 NaOH vs uninjured ===\n"); print(head(net, 20))

## ---------- 1b. 细胞类型层面的总通讯强度 ----------
ct_strength <- function(cc, tag) {
  P <- cc@netP$prob
  tot <- apply(P, c(1, 2), sum, na.rm = TRUE)   # source x target
  dt <- as.data.table(as.table(tot)); setnames(dt, c("source", "target", "strength"))
  dt[, run := tag][]
}
cs <- rbind(ct_strength(cc_naoh, "NaOH"), ct_strength(cc_un, "uninjured"))
cs_w <- dcast(cs, source + target ~ run, value.var = "strength", fill = 0)
setnames(cs_w, c("NaOH", "uninjured"), c("strength_NaOH", "strength_uninjured"))
fwrite(cs_w, file.path(RES, "D9_cellchat_celltype_strength.tsv"), sep = "\t")
cat("\n=== 细胞对通讯强度 (NaOH) top12 ===\n"); print(head(cs_w[order(-strength_NaOH)], 12))

## ---------- 2. 炎症细胞 → 血管细胞 ----------
inflam <- c("Neutrophil", "Macrophage")
vessel <- c("Endothelium", "Endothelium_Lymphatic")
extract_st <- function(cc, tag) {
  P <- cc@netP$prob
  rbindlist(lapply(dimnames(P)[[3]], function(k) {
    m <- P[, , k]
    dt <- as.data.table(as.table(m)); setnames(dt, c("source", "target", "prob"))
    dt[, pathway := k][source %in% inflam & target %in% vessel]
  }))[, run := tag]
}
st <- rbind(extract_st(cc_naoh, "NaOH"), extract_st(cc_un, "uninjured"))
agg <- st[, .(prob = sum(prob, na.rm = TRUE)), by = .(run, source, target, pathway)][order(-prob)]
fwrite(agg, file.path(RES, "D9_cellchat_inflam_to_vessel_real.tsv"), sep = "\t")
cat("\n=== 炎症 → 血管（NaOH）top15 ===\n")
print(head(agg[run == "NaOH"], 15))

## ---------- 3. LR 对层面（net$pval 的第 3 维是 LR 对字符串） ----------
sig <- cc_naoh@net$pval
sig_pos <- which(sig < 0.05 & sig > -1, arr.ind = TRUE)
lrdt <- data.table(source = dimnames(sig)[[1]][sig_pos[, 1]],
                   target = dimnames(sig)[[2]][sig_pos[, 2]],
                   lr     = dimnames(sig)[[3]][sig_pos[, 3]],
                   pval   = sig[sig_pos])
lrdt <- lrdt[source %in% inflam & target %in% vessel][order(pval)]
prob_long <- data.table(source = dimnames(cc_naoh@net$prob)[[1]][sig_pos[, 1]],
                        target = dimnames(cc_naoh@net$prob)[[2]][sig_pos[, 2]],
                        lr     = dimnames(cc_naoh@net$prob)[[3]][sig_pos[, 3]],
                        prob   = cc_naoh@net$prob[sig_pos])
lrdt <- merge(lrdt, prob_long, by = c("source", "target", "lr"))[order(-prob)]
fwrite(lrdt, file.path(RES, "D9_cellchat_lr_to_vessel_real.tsv"), sep = "\t")
cat("\n=== 炎症 → 血管 显著 LR（p<0.05）条目数:", nrow(lrdt), "===\n")
if (nrow(lrdt)) print(head(lrdt, 18))

## ---------- 4. 图 ----------
# (a) 通路气泡图
pdat <- head(net[strength_NaOH > 0][order(-strength_NaOH)], 22)
pdat[, pathway := factor(pathway, levels = rev(pathway))]
pa <- ggplot(pdat, aes(x = log2FC, y = pathway, size = strength_NaOH)) +
  geom_vline(xintercept = 0, linewidth = 0.3, linetype = 2) +
  geom_point(aes(colour = log2FC), alpha = 0.85) +
  scale_colour_gradient2(low = "#2166AC", mid = "grey70", high = "#B2182B") +
  scale_size_continuous(range = c(0.6, 4.2)) +
  labs(x = expression(log[2]~"(NaOH / uninjured)"), y = NULL, size = "strength",
       title = "CellChat pathway comparison (NaOH vs uninjured)",
       subtitle = "Mouse CellChatDB, trimean model, permutation-tested") +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = 6.6, colour = "grey30"))
ggsave(file.path(FIG, "D9_Fig9a_cellchat_pathways.png"), pa, width = 5.2, height = 4.4, dpi = 300, bg = "white")
ggsave(file.path(FIG, "D9_Fig9a_cellchat_pathways.pdf"), pa, width = 5.2, height = 4.4)

# (b) 炎症 → 血管 分离条形（source x target 分面）
pv <- agg[run == "NaOH"][order(-prob)]
pv <- head(pv, 16)
pv[, pair := paste(source, "→", target)]
pv[, pair := factor(pair, levels = rev(unique(pair)))]
pb <- ggplot(pv, aes(x = prob, y = pair, fill = pathway)) +
  geom_col(width = 0.72, linewidth = 0) +
  labs(x = "Communication probability", y = NULL, fill = "pathway",
       title = "Inflammatory to vascular signalling, alkali burn",
       subtitle = "Real CellChat probabilities, top pathway-cell pairs") +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        legend.key.size = unit(0.32, "cm"),
        legend.text = element_text(size = 6.2),
        plot.subtitle = element_text(size = 6.6, colour = "grey30"))
ggsave(file.path(FIG, "D9_Fig9b_cellchat_inflam_vessel.png"), pb, width = 5.8, height = 3.9, dpi = 300, bg = "white")
ggsave(file.path(FIG, "D9_Fig9b_cellchat_inflam_vessel.pdf"), pb, width = 5.8, height = 3.9)

cat("\n=== DONE 06 ===\n")