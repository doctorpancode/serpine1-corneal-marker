# ============================================================
# D4-05 RFE 简历汇总 + 重绘 Fig5（列名: Normal/Lesion/Overall/var/Variables/Resample）
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
out <- "D:/A生信分析/虚拟敲除/第4天执行结果"
sv <- fread(file.path(out, "results", "D4_svmrfe_profile.tsv"))
cat("columns:", paste(names(sv), collapse = ", "), "\n")

# 每个 (Variables, Resample) 的整体准确率 = 该折该子集内 Overall 的均值
agg <- sv[, .(acc = mean(Overall, na.rm = TRUE)), by = .(Variables, Resample)]
pv <- agg[, .(Accuracy = mean(acc), AccuracySD = sd(acc)), by = Variables][order(Variables)]
fwrite(pv, file.path(out, "results", "D4_svmrfe_accuracy.tsv"), sep = "\t")
print(pv)

best_n <- 8
p3 <- ggplot(pv, aes(Variables, Accuracy)) +
  geom_errorbar(aes(ymin = Accuracy - AccuracySD, ymax = Accuracy + AccuracySD),
                width = 0.03, colour = "grey60") +
  geom_line(colour = "#4c72b0") + geom_point(colour = "#4c72b0", size = 2.2) +
  geom_point(data = pv[Variables == best_n], colour = "#d62728", size = 3.8) +
  scale_x_continuous(trans = "log2", breaks = pv$Variables) +
  labs(title = "SVM-RFE recursive feature elimination (5-fold CV)",
       subtitle = paste0("red = ", best_n, "-gene minimal subset reaching CV accuracy ",
                         sprintf("%.2f", pv[Variables == best_n]$Accuracy),
                         " (n = 21; upper-bound plateau)"),
       x = "Number of retained features (log2)", y = "Mean CV accuracy") +
  theme_bw(base_size = 11)
ggsave(file.path(out, "figures", "D4_Fig5_svmrfe.pdf"), p3, width = 6.6, height = 4.7)
ggsave(file.path(out, "figures", "D4_Fig5_svmrfe.png"), p3, width = 6.6, height = 4.7, dpi = 300)
cat("Fig5 redrawn\n")