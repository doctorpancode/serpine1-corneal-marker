# ============================================================
# D9-04  待办 10.3-(2)：多敲除（multi-KO）
#
#   官方 hkManifoldAlignment 只支持"每次敲一个"。但 .kernelEffect 的形式是
#       effect(gKO) = sum_{g in gKO}  w_g * H[g, ]  * sd
#     w_g = -mu_g / sd_g
#   故对多个基因同时敲除，扰动向量 = 各单敲除"核效应向量"的**线性叠加**
#   （因为 H 是全矩阵，行与行可以直接相加；著名的非线性只来自后续 boxcox+Z 标准化）。
#   本脚本：
#    (a) 直接叠加：D_multi = sum of |kernel vectors|（保留符号场再取绝对值）
#    (b) 交互效应：D_multi - (D_a + D_b)，若显著超出即存在"协同/拮抗"
#    (c) 检验 SERPINE1 + TMSB4X 是否真的"相反两端"：
#        sign(kernel vector) 的一致性（同向率）
#    (d) 三敲除组合（Serpine1 + Tmsb4x + P4hb）
#
#   产出：
#     results/D9_multiko_scores.tsv
#     results/D9_multiko_summary.tsv
#     figures/main/D9_Fig7_multiko.{png,pdf}
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(scTenifoldKnk)
})
OUT <- "D:/A生信分析/虚拟敲除/第9天执行结果"
D8  <- "D:/A生信分析/虚拟敲除/第8天执行结果"
RES <- file.path(OUT, "results"); FIG <- file.path(OUT, "figures", "main")

CORE <- c("Serpine1","Tmsb4x","P4hb","S100a13","Vim")
EPI  <- c("Krt12","Aqp5","Gsta4","Upk3bl","Upk1b","Krt15","Krt17","Krt19","Krt13",
          "Dsg1a","Perp","Tmprss11d","Lypd3","Calml3","Gsta2","Sfn","Krt5","Krt14")

WT  <- as.matrix(readRDS(file.path(D8, "data", "D8_WT_network.rds")))
qcX <- readRDS(file.path(D8, "data", "D8_knk_qcCounts.rds"))
diag(WT) <- 0
genes <- rownames(WT)
H <- scTenifoldKnk::heatKernel(WT, t = 10)
S <- scTenifoldKnk:::.logCPM(as(qcX, "dgCMatrix"), genes)
EPI2 <- intersect(EPI, genes)

# 单个 (或一组) 基因的核效应向量（带符号）
kernel_vec <- function(gk) {
  w <- -S$mu[gk] / S$sd[gk]
  as.numeric(crossprod(w, H[gk, , drop = FALSE])) * S$sd
}
names_kernel <- function(gk) { v <- kernel_vec(gk); names(v) <- genes; v }

combos <- list(
  SERPINE1 = "Serpine1",
  TMSB4X   = "Tmsb4x",
  P4HB     = "P4hb",
  `SERPINE1+TMSB4X` = c("Serpine1","Tmsb4x"),
  `SERPINE1+P4HB`   = c("Serpine1","P4hb"),
  `TMSB4X+P4HB`     = c("Tmsb4x","P4hb"),
  `SERPINE1+TMSB4X+P4HB` = c("Serpine1","Tmsb4x","P4hb")
)

## ---------- 1. 叠加 vs 加性 ----------
rows <- list()
for (nm in names(combos)) {
  gk <- combos[[nm]]
  v_multi <- kernel_vec(gk)                    # 叠加（线性）
  v_add <- rowSums(sapply(gk, function(g) kernel_vec(g)))  # 与叠加应完全相同（数值校验）
  stopifnot(max(abs(v_multi - v_add)) < 1e-10 * max(abs(v_multi)))
  # 与"逐基因单独敲除后取绝对值再取平均"作对比（非线性的直觉版本）
  v_indiv <- rowMeans(abs(sapply(gk, function(g) kernel_vec(g))))
  epi_multi <- mean(abs(v_multi[EPI2]))
  epi_indiv <- mean(v_indiv[EPI2])
  rows[[nm]] <- data.table(combo = nm, n_genes = length(gk),
                           epi_mean_abs_multi = epi_multi,
                           epi_mean_abs_individual_avg = epi_indiv,
                           ratio = epi_multi / epi_indiv,
                           top_gene = genes[which.max(abs(v_multi))],
                           top_abs = max(abs(v_multi)))
}
mk <- rbindlist(rows)
fwrite(mk, file.path(RES, "D9_multiko_summary.tsv"), sep = "\t")
cat("\n=== 多敲除组合 ===\n"); print(mk)

## ---------- 2. SERPINE1 vs TMSB4X 的方向一致性 ----------
v_sp <- names_kernel("Serpine1"); v_tx <- names_kernel("Tmsb4x")
same_sign <- mean(sign(v_sp[EPI2]) == sign(v_tx[EPI2]))
cat("\n上皮 18 基因中 SERPINE1 与 TMSB4X 核效应同号比例:", round(same_sign, 3), "\n")
ct <- cor.test(v_sp, v_tx, method = "pearson")
cat("全 947 基因核效应向量 Pearson r =", round(ct$estimate, 4), " p =", signif(ct$p.value, 3), "\n")
dir_dt <- data.table(gene = genes,
                     kernel_SERPINE1 = v_sp,
                     kernel_TMSB4X = v_tx,
                     on_epi = genes %in% EPI2)
fwrite(dir_dt, file.path(RES, "D9_multiko_kernel_vectors.tsv"), sep = "\t")

## ---------- 3. 交互项（叠加 − 各自绝对值之和） ----------
inter <- abs(v_sp + v_tx) - (abs(v_sp) + abs(v_tx))   # <=0 恒成立（三角不等式）
inter_ratio <- mean(abs(v_sp[EPI2]) + abs(v_tx[EPI2])) / mean(abs(v_sp[EPI2] + v_tx[EPI2]))
cat("上皮区间 |sp|+|tx| / |sp+tx| =", round(inter_ratio, 3),
    "（>1 说明两者在上皮上部分抵消）\n")

summ2 <- data.table(
  metric = c("epi_same_sign_fraction", "pearson_r_all_genes",
             "p_all_genes", "cancellation_index_epi"),
  value = c(same_sign, unname(ct$estimate), ct$p.value, inter_ratio)
)
fwrite(summ2, file.path(RES, "D9_multiko_direction.tsv"), sep = "\t")
print(summ2)

## ---------- 4. 图 D9-Fig7：单敲 vs 叠加（上皮基因条形对照） ----------
bar <- data.table(gene = EPI2,
                  SERPINE1 = v_sp[EPI2],
                  TMSB4X   = v_tx[EPI2],
                  `SERPINE1+TMSB4X` = (v_sp + v_tx)[EPI2])
bl <- melt(bar, id.vars = "gene", variable.name = "knockout", value.name = "kernel_effect")
bl[, knockout := factor(knockout, levels = c("SERPINE1","TMSB4X","SERPINE1+TMSB4X"))]
bl[, gene := factor(gene, levels = bar[order(SERPINE1)]$gene)]
p7 <- ggplot(bl, aes(x = gene, y = kernel_effect, fill = knockout)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.72, linewidth = 0.1) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c(SERPINE1 = "#B2182B", TMSB4X = "#2166AC",
                               `SERPINE1+TMSB4X` = "#762A83")) +
  labs(x = NULL, y = "Kernel effect on gene (signed)",
       title = "Single vs combined in-silico knockout of SERPINE1 and TMSB4X",
       subtitle = sprintf("Effects add linearly in the heat-kernel; on the epithelial module the two vectors partially cancel (index %.2f).  Same-sign fraction = %.2f",
                          inter_ratio, same_sign),
       fill = NULL) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        panel.grid.minor = element_blank(),
        legend.position = "top",
        plot.subtitle = element_text(size = 6.8, colour = "grey30"))
ggsave(file.path(FIG, "D9_Fig7_multiko.png"), p7, width = 7.4, height = 3.8, dpi = 300)
ggsave(file.path(FIG, "D9_Fig7_multiko.pdf"), p7, width = 7.4, height = 3.8)
cat("\n[saved] D9_Fig7_multiko\n")

## ---------- 5. 散点：两向量在全基因层面的关系 ----------
p7b <- ggplot(dir_dt, aes(x = kernel_SERPINE1, y = kernel_TMSB4X)) +
  geom_point(aes(colour = on_epi), size = 0.75, alpha = 0.6) +
  geom_hline(yintercept = 0, linewidth = 0.25, linetype = 2) +
  geom_vline(xintercept = 0, linewidth = 0.25, linetype = 2) +
  scale_colour_manual(values = c(`FALSE` = "grey65", `TRUE` = "#D95F02"),
                      labels = c("other genes", "epithelial module")) +
  labs(x = "Kernel effect of SERPINE1 KO", y = "Kernel effect of TMSB4X KO",
       colour = NULL,
       title = "Opposite-sign structure of the two knockouts",
       subtitle = sprintf("Pearson r = %.3f across all %d genes", unname(ct$estimate), length(genes))) +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        plot.subtitle = element_text(size = 7, colour = "grey30"))
ggsave(file.path(FIG, "D9_Fig7b_multiko_scatter.png"), p7b, width = 4.0, height = 3.8, dpi = 300)
ggsave(file.path(FIG, "D9_Fig7b_multiko_scatter.pdf"), p7b, width = 4.0, height = 3.8)
cat("[saved] D9_Fig7b_multiko_scatter\n")

cat("\n=== DONE 04 ===\n")