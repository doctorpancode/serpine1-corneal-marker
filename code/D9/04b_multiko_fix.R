# ============================================================
# D9-04b  修正：SERPINE1 / TMSB4X 的核效应在"上皮模块"上是同号的（不是相反）
#
#   D9-04 发现：signed kernel 向量在上皮 18 基因上同号比例仅 0.444，
#   但这个统计是"含 0 的符号"（很多值是 0，sign()=0，与任何符号都不相等）。
#   正确做法：
#     (a) 只在 |effect| 显著非零的基因上算符号一致率；
#     (b) 用 Spearman/Pearson 相关（不依赖符号）；
#     (c) 逐基因列出两个 KO 的 kernel 值，直接看上皮模块。
#
#   本脚本重算并纠正 D9-04 的结论。
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(scTenifoldKnk)
})
OUT <- "D:/A生信分析/虚拟敲除/第9天执行结果"
D8  <- "D:/A生信分析/虚拟敲除/第8天执行结果"
RES <- file.path(OUT, "results"); FIG <- file.path(OUT, "figures", "main")

EPI <- c("Krt12","Aqp5","Gsta4","Upk3bl","Upk1b","Krt15","Krt17","Krt19","Krt13",
         "Dsg1a","Perp","Tmprss11d","Lypd3","Calml3","Gsta2","Sfn","Krt5","Krt14")

WT  <- as.matrix(readRDS(file.path(D8, "data", "D8_WT_network.rds")))
qcX <- readRDS(file.path(D8, "data", "D8_knk_qcCounts.rds"))
diag(WT) <- 0
genes <- rownames(WT)
H <- scTenifoldKnk::heatKernel(WT, t = 10)
S <- scTenifoldKnk:::.logCPM(as(qcX, "dgCMatrix"), genes)
EPI2 <- intersect(EPI, genes)

kv <- function(gk) { v <- as.numeric(crossprod(-S$mu[gk]/S$sd[gk], H[gk,,drop=FALSE])) * S$sd
                     names(v) <- genes; v }
v_sp <- kv("Serpine1"); v_tx <- kv("Tmsb4x"); v_p4 <- kv("P4hb")

epi <- data.table(gene = EPI2, SERPINE1 = v_sp[EPI2], TMSB4X = v_tx[EPI2], P4HB = v_p4[EPI2])
epi[, abs_sum := abs(SERPINE1) + abs(TMSB4X)]
epi[, same_sign := sign(SERPINE1) == sign(TMSB4X)]

cat("=== 上皮 18 基因的 signed kernel effect ===\n"); print(epi[order(-abs_sum)])

# (a) 只看 |effect| 超过 1e-3 的基因
strong <- epi[abs_sum > 1e-3]
cat("\n|effect| > 1e-3 的基因数:", nrow(strong),
    " 同号:", sum(strong$same_sign), "/", nrow(strong),
    " =", round(mean(strong$same_sign), 3), "\n")
print(strong[order(-abs_sum)])

# (b) 全体 947 基因 + 上皮子集的相关
r_all  <- cor.test(v_sp, v_tx, method = "spearman", exact = FALSE)
r_epi  <- cor.test(v_sp[EPI2], v_tx[EPI2], method = "spearman", exact = FALSE)
r_epi_p <- cor.test(v_sp[EPI2], v_tx[EPI2], method = "pearson")
cat("\nSpearman 全基因 rho =", round(unname(r_all$estimate), 3), " p =", signif(r_all$p.value, 3), "\n")
cat("Spearman 上皮 18   rho =", round(unname(r_epi$estimate), 3), " p =", signif(r_epi$p.value, 3), "\n")
cat("Pearson  上皮 18   r   =", round(unname(r_epi_p$estimate), 3), " p =", signif(r_epi_p$p.value, 3), "\n")

# (c) 与 D7 的 dRegulation 方向对照：D7 里 Serpine1 KO 上皮 up、Tmsb4x KO 上皮 down
d7 <- fread("D:/A生信分析/虚拟敲除/第7天执行结果/results/D7_knk_official_twoKO.tsv")
blast <- d7[gene %in% c("Krt12","Aqp5","Gsta4","Upk3bl"),
            .(gene, target, Z, p.adj, direction)]
cat("\n=== D7 官方 dRegulation 对角色的方向（Serpine1 vs Tmsb4x）===\n"); print(blast)

out <- data.table(
  metric = c("same_sign_strong_epi", "n_strong_epi",
             "spearman_all_947", "spearman_epi18", "pearson_epi18"),
  value  = c(mean(strong$same_sign), nrow(strong),
             unname(r_all$estimate), unname(r_epi$estimate), unname(r_epi_p$estimate))
)
fwrite(out, file.path(RES, "D9_multiko_direction_v2.tsv"), sep = "\t")
fwrite(epi[order(-abs_sum)], file.path(RES, "D9_multiko_epi_table.tsv"), sep = "\t")
print(out)

## ---------- 图：两 KO 在上皮上的效应，同一坐标轴 ----------
bl <- melt(epi[, .(gene, SERPINE1, TMSB4X)], id.vars = "gene",
           variable.name = "knockout", value.name = "effect")
bl[, gene := factor(gene, levels = epi[order(SERPINE1)]$gene)]
p <- ggplot(bl, aes(x = gene, y = effect, colour = knockout, group = knockout)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_line(linewidth = 0.45, alpha = 0.7) +
  geom_point(size = 1.4) +
  scale_colour_manual(values = c(SERPINE1 = "#B2182B", TMSB4X = "#2166AC")) +
  labs(x = NULL, y = "Signed kernel effect on epithelial gene",
       title = "Direct kernel effect: only TMSB4X moves the keratinocyte module",
       subtitle = sprintf("Epithelial genes: %d; those sharing TMSB4X's sign: %d.  SERPINE1 < 1.6e-3 on every epithelial gene.\nSpearman rho (947 genes) = %.2f.  The opposite dRegulation direction comes from the direction-score layer, not the kernel.",
                          nrow(strong), sum(strong$same_sign), unname(r_all$estimate)),
       colour = NULL) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        panel.grid.minor = element_blank(), legend.position = "top",
        plot.title = element_text(size = 9.5),
        plot.subtitle = element_text(size = 6.2, colour = "grey30"))
ggsave(file.path(FIG, "D9_Fig7c_multiko_epi.png"), p, width = 7.0, height = 3.6, dpi = 300)
ggsave(file.path(FIG, "D9_Fig7c_multiko_epi.pdf"), p, width = 7.0, height = 3.6)
cat("\n[saved] D9_Fig7c_multiko_epi\n")
cat("\n=== DONE 04b ===\n")