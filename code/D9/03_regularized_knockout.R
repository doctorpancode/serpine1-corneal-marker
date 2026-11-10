# ============================================================
# D9-03  待办 10.3-(1)：正则化虚拟敲除 —— 把中心性从扰动谱里剥掉
#
#   问题（D8 已知）：947 个 KO 的上皮 mean|Δ| 与网络出度 Spearman ρ=0.9169。
#   说明 hkManifoldAlignment 的 kernelEffect 本质是
#       w = -mu[ko]/sd[ko];  effect = crossprod(w, H[ko,]) * sd
#   即"敲除基因越中心（H 行越重），谱越宽"——这是拓扑产物，不是生物学。
#
#   本脚本实现三种正则化：
#    (R1) 行归一：把 H 的每一行除以该行 L1 范数（把"权重总量"归一化）
#    (R2) 谱归一：把每个 KO 的扰动向量按其 2-范数归一后取方向（去掉幅度）
#    (R3) 中心性残差：用 loess(mean|Δ| ~ out_degree) 取残差，得到
#         "超出其拓扑预期的上皮扰动"
#
#   评价标准（这是 10.3 里最关键的判据）：
#    - 正则化后 ρ(out_degree, mean|Δ|) 必须显著下降；
#    - TMSB4X 仍居前（说明它是真信号），SERPINE1 排名若上升，
#      则"中心性混杂"被部分剥离，可在讨论里写成"正则化后 SERPINE1 上升 至第 N 位"。
#
#   产出：
#     results/D9_regularized_ko_scores.tsv
#     results/D9_regularization_compare.tsv
#     figures/main/D9_Fig6_regularization.{png,pdf}
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(scTenifoldKnk)
})
OUT <- "D:/A生信分析/虚拟敲除/第9天执行结果"
D8  <- "D:/A生信分析/虚拟敲除/第8天执行结果"
RES <- file.path(OUT, "results"); FIG <- file.path(OUT, "figures", "main")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)

CORE <- c("Serpine1","Tmsb4x","P4hb","S100a13","Sod1","Prdx2","Vdac1","Capg","Gja1","Xiap","Vim")
EPI  <- c("Krt12","Aqp5","Gsta4","Upk3bl","Upk1b","Krt15","Krt17","Krt19","Krt13",
          "Dsg1a","Perp","Tmprss11d","Lypd3","Calml3","Gsta2","Sfn","Krt5","Krt14")

## ---------- 1. 重建 H（热核）与表达统计 ----------
WT  <- as.matrix(readRDS(file.path(D8, "data", "D8_WT_network.rds")))
qcX <- readRDS(file.path(D8, "data", "D8_knk_qcCounts.rds"))
diag(WT) <- 0
genes <- rownames(WT)
cat("WT:", dim(WT), "\n")

H <- scTenifoldKnk::heatKernel(WT, t = 10)
cat("H:", dim(H), " 行和范围:", paste(round(range(rowSums(H)), 4), collapse = " ~ "), "\n")

# .logCPM 统计（与官方 .kernelEffect 同源的 mu/sd）
S <- scTenifoldKnk:::.logCPM(as(qcX, "dgCMatrix"), genes)
cat("S$mu 长度:", length(S$mu), "\n")

EPI2 <- intersect(EPI, genes)
cat("参与评价的上皮基因:", length(EPI2), "/", length(EPI), "\n")

## ---------- 2. 三种正则化 ----------
# 原始效应（复现官方）：effect_i = (-mu_i/sd_i) %*% H[i,] * sd
raw_eff <- function(i) abs(as.numeric(crossprod(-S$mu[i] / S$sd[i], H[i, ])) * S$sd)
# R1：H 行 L1 归一
H1 <- H / rowSums(abs(H))
r1_eff <- function(i) abs(as.numeric(crossprod(-S$mu[i] / S$sd[i], H1[i, ])) * S$sd)

n <- length(genes)
E_raw <- matrix(NA_real_, n, n, dimnames = list(genes, genes))
E_r1  <- matrix(NA_real_, n, n, dimnames = list(genes, genes))
pb <- txtProgressBar(max = n, style = 3)
for (i in seq_len(n)) {
  E_raw[i, ] <- raw_eff(i)
  E_r1[i, ]  <- r1_eff(i)
  setTxtProgressBar(pb, i)
}
close(pb)
diag(E_raw)[] <- NA; diag(E_r1)[] <- NA
cat("\n原始/行归一阵列计算完成\n")

outdeg <- rowSums(WT != 0)
base <- data.table(ko_gene = genes, out_degree = as.numeric(outdeg))
base[, mean_raw   := rowMeans(E_raw[, EPI2], na.rm = TRUE)]
base[, mean_r1    := rowMeans(E_r1[,  EPI2], na.rm = TRUE)]
# R2：谱方向余弦（把每个 KO 的上皮剖面归一为方向）
epi_mat <- E_raw[, EPI2]
norm2 <- sqrt(rowSums(epi_mat^2, na.rm = TRUE))
dir_cos <- epi_mat / norm2                       # 单位向量
# 以 TMSB4X 的方向为参考，算与它的余弦相似度（生物学意义上的"同向扰动"）
ref <- dir_cos["Tmsb4x", ]
base[, cos_to_tmsb4x := as.numeric(dir_cos %*% ref)]
# R3：中心性残差
fit_lo <- loess(mean_raw ~ out_degree, data = base, span = 0.5)
base[, resid_raw := resid(fit_lo)]
fit_lo2 <- loess(mean_r1 ~ out_degree, data = base, span = 0.5)
base[, resid_r1 := resid(fit_lo2)]

## ---------- 3. 评价：ρ 是否被压下去 ----------
sp <- function(x, y) suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
r0 <- sp(base$out_degree, base$mean_raw)
r1 <- sp(base$out_degree, base$mean_r1)
cat("\n原始 ρ(out_degree, mean_raw) =", round(unname(r0$estimate), 4),
    "  p =", signif(r0$p.value, 3), "\n")
cat("行归一 ρ(out_degree, mean_r1) =", round(unname(r1$estimate), 4),
    "  p =", signif(r1$p.value, 3), "\n")

# 排名
base[, rank_raw  := frank(-mean_raw,  ties.method = "first")]
base[, rank_r1   := frank(-mean_r1,   ties.method = "first")]
base[, rank_cos  := frank(-cos_to_tmsb4x, ties.method = "first")]
base[, rank_resid := frank(-resid_raw, ties.method = "first")]

cmp <- base[ko_gene %in% CORE]
setorder(cmp, rank_raw)
fwrite(base, file.path(RES, "D9_regularized_ko_scores.tsv"), sep = "\t")

summ <- data.table(
  metric = c("rho_outdegree_vs_mean_raw",
             "rho_outdegree_vs_mean_rowL1norm",
             "rho_outdegree_vs_cos_to_TMSB4X"),
  rho = c(unname(r0$estimate), unname(r1$estimate),
          unname(sp(base$out_degree, base$cos_to_tmsb4x)$estimate)),
  p = c(r0$p.value, r1$p.value, sp(base$out_degree, base$cos_to_tmsb4x)$p.value)
)
fwrite(summ, file.path(RES, "D9_regularization_compare.tsv"), sep = "\t")
cat("\n=== 三种度量的中心性相关 ===\n"); print(summ)

cat("\n=== 核心基因在四种度量下的排名（947 KO）===\n")
print(cmp[, .(ko_gene, out_degree, mean_raw = signif(mean_raw, 4),
              rank_raw, rank_r1, cos_to_tmsb4x = round(cos_to_tmsb4x, 3), rank_cos,
              resid_raw = signif(resid_raw, 4), rank_resid)])

## ---------- 4. 图 D9-Fig6：正则化前后对比 ----------
pd <- rbind(
  base[, .(ko_gene, out_degree, value = mean_raw, metric = "Raw |delta| (official)")],
  base[, .(ko_gene, out_degree, value = mean_r1,  metric = "Row-L1-normalised |delta|")]
)
pd[, metric := factor(metric, levels = c("Raw |delta| (official)", "Row-L1-normalised |delta|"))]
pd[, is_core := ko_gene %in% CORE]
core_pd <- pd[is_core == TRUE]

p6 <- ggplot(pd, aes(x = out_degree, y = value)) +
  geom_point(size = 0.45, alpha = 0.25, colour = "grey40") +
  geom_smooth(method = "loess", span = 0.6, se = FALSE, linewidth = 0.6, colour = "#2166AC") +
  geom_point(data = core_pd, size = 1.9, shape = 21, fill = "#B2182B",
             colour = "black", stroke = 0.25) +
  ggrepel::geom_text_repel(data = core_pd, aes(label = ko_gene), size = 2.1,
                           box.padding = 0.28, min.segment.length = 0,
                           segment.size = 0.18, max.overlaps = 30) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(x = "Network out-degree", y = "Epithelial mean |delta| (18 genes)",
       title = "Regularisation does not remove the centrality bias",
       subtitle = sprintf("Spearman rho = %.3f (raw)  ->  %.3f (row-L1-normalised);  n = %d knockouts",
                          unname(r0$estimate), unname(r1$estimate), nrow(base))) +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "none",
        strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text = element_text(face = "bold", size = 8),
        plot.subtitle = element_text(size = 7, colour = "grey30"))
ggsave(file.path(FIG, "D9_Fig6_regularization.png"), p6, width = 8.4, height = 3.6, dpi = 300)
ggsave(file.path(FIG, "D9_Fig6_regularization.pdf"), p6, width = 8.4, height = 3.6)
cat("\n[saved] D9_Fig6_regularization\n")
cat("\n=== DONE 03 ===\n")