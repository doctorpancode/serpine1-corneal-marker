# 用 R 的 rsvg 把 SVG 栅格化成 PNG 预览
out <- "D:/A生信分析/虚拟敲除/第8天执行结果/figures/editable"
src <- file.path(out, "Fig9_mechanism_model_editable.svg")
rsvg::rsvg_png(src, file.path(out, "Fig9_mechanism_model_editable.png"), width = 1520)
cat("png written\n")

