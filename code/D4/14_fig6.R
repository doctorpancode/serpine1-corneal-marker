# ============================================================
# D4-14 Fig6 三算法证据矩阵（04_converge.R 因列名错误中断，此脚本单独补出）
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
out <- "D:/A生信分析/虚拟敲除/第4天执行结果"
vote <- fread(file.path(out, "results", "D4_ml_vote_matrix.tsv"))
sel <- vote[n_vote > 0][order(-n_vote, -PPI_hub, -rf_MDG)]
sel <- head(sel, 30)
vm <- melt(sel, id.vars = "gene", measure.vars = c("LASSO", "RF_top20", "SVM_RFE"),
           variable.name = "method", value.name = "hit")
vm[, gene := factor(gene, levels = rev(sel$gene))]
vm[, method := factor(method, levels = c("LASSO", "RF_top20", "SVM_RFE"),
                      labels = c("LASSO", "RF top20", "SVM-RFE"))]
lab <- sel[, .(gene = factor(gene, levels = rev(sel$gene)), n = n_vote, hub = PPI_hub)]
p <- ggplot(vm, aes(method, gene, fill = factor(hit))) +
  geom_tile(colour = "white", linewidth = 0.9) +
  scale_fill_manual(values = c("0" = "grey93", "1" = "#d62728"), guide = "none") +
  geom_text(data = lab, aes(x = 3.62, y = gene, label = n), inherit.aes = FALSE,
            size = 3, fontface = "bold") +
  labs(title = "D4 three-algorithm evidence matrix (top 30 by vote count)",
       subtitle = "right column = number of algorithms selecting the gene; only TMSB4X reaches 3/3",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size = 7.5))
ggsave(file.path(out, "figures", "D4_Fig6_ml_convergence.pdf"), p, width = 6.4, height = 6.4)
ggsave(file.path(out, "figures", "D4_Fig6_ml_convergence.png"), p, width = 6.4, height = 6.4, dpi = 300)
cat("Fig6 done\n")