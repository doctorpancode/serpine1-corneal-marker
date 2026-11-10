# ============================================================
# D4-11 决策矩阵：把 PPI / 三算法 / 跨队列 / 生物学 四维证据汇总成一张表
# 输出: results/D4_candidate_scorecard.tsv
#       figures/D4_Fig9_scorecard.pdf/.png
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
out <- "D:/A生信分析/虚拟敲除/第4天执行结果"

vote <- fread(file.path(out, "results", "D4_ml_vote_matrix.tsv"))
hub  <- fread(file.path(out, "results", "D4_hub_intersection.tsv"))
r2   <- fread(file.path(out, "results", "D4_crosscohort_final.tsv"))
w3   <- fread(file.path(out, "results", "D4_crosscohort_mouse_final.tsv"))

cand <- c("TMSB4X", "SERPINE1", "S100A13", "XIAP", "CAPG", "PRDX2", "VDAC1", "P4HB", "SOD1", "GJA1")
sc <- data.table(gene = cand)
sc <- merge(sc, vote[, .(gene, LASSO, RF_top20, SVM_RFE, PPI_hub, n_vote, lasso_coef)],
            by = "gene", all.x = TRUE)

# 人角膜跨队列：取"病变方向最一致"的显著分层
h <- r2[gene %in% cand, .(human_best_p = min(p),
                          human_sig = any(p < 0.05),
                          human_direction = paste(unique(fifelse(delta > 0, "up", "down")), collapse = "/")),
        by = gene]
sc <- merge(sc, h, by = "gene", all.x = TRUE)
m <- w3[gene %in% cand, .(mouse_d5 = log2FC_d5, mouse_d15 = log2FC_d15), by = gene]
sc <- merge(sc, m, by = "gene", all.x = TRUE)
sc[, mouse_peak_direction := fifelse(abs(mouse_d5) >= abs(mouse_d15),
                                     fifelse(mouse_d5 > 0, "up", "down"),
                                     fifelse(mouse_d15 > 0, "up", "down"))]
sc[, mouse_magnitude := pmax(abs(mouse_d5), abs(mouse_d15))]
sc[, human_dir_clean := fifelse(!is.na(human_direction) & human_direction %in% c("up", "down"),
                                human_direction, NA_character_)]
sc[, concordant := human_dir_clean == mouse_peak_direction]

# 综合分（0-8）：三算法票(0-3) + PPI hub(0-2) + 人角膜 p<0.05(0-1) + 鼠时序大幅变化(>1.5 log2 -> 0-2)
sc[, score := n_vote + 2 * PPI_hub + as.integer(human_sig) +
     2 * as.integer(mouse_magnitude > 1.5)]
setorder(sc, -score, -n_vote, human_best_p)

sc[, prior := gene %in% c("SERPINE1", "THBS1", "LOX")]
fwrite(sc, file.path(out, "results", "D4_candidate_scorecard.tsv"), sep = "\t")
print(sc[, .(gene, n_vote, PPI_hub, human_best_p = signif(human_best_p, 3), human_sig,
             human_dir_clean, mouse_peak_direction, mouse_magnitude = round(mouse_magnitude, 2),
             concordant, score)])

# 热图式 scorecard
m1 <- melt(sc, id.vars = "gene",
           measure.vars = c("LASSO", "RF_top20", "SVM_RFE", "PPI_hub", "human_sig"),
           variable.name = "evidence", value.name = "present")
m1[, gene := factor(gene, levels = rev(sc$gene))]
p <- ggplot(m1, aes(evidence, gene, fill = factor(present))) +
  geom_tile(colour = "white", linewidth = 0.8) +
  scale_fill_manual(values = c("0" = "grey93", "1" = "#d62728"),
                    labels = c("absent", "present"), name = NULL) +
  geom_text(data = sc, aes(x = 5.75, y = gene, label = sprintf("%.1f", score)),
            inherit.aes = FALSE, size = 3.2, fontface = "bold") +
  labs(title = "D4 candidate scorecard: four independent evidence layers",
       subtitle = "columns = LASSO / RF / SVM-RFE / cytoHubba intersection / human-cornea significance; right = composite score",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid = element_blank())
ggsave(file.path(out, "figures", "D4_Fig9_scorecard.pdf"), p, width = 8.2, height = 4.8)
ggsave(file.path(out, "figures", "D4_Fig9_scorecard.png"), p, width = 8.2, height = 4.8, dpi = 300)
cat("\nD4-11 done\n")