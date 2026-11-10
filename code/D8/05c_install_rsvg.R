# 尝试安装 rsvg（CRAN 有 Windows 二进制包），用于把可编辑 SVG 栅格化预览
options(repos = c(CRAN = "https://cloud.r-project.org"))
tryCatch({
  install.packages("rsvg", dependencies = FALSE, quiet = TRUE)
  cat("rsvg installed:", requireNamespace("rsvg", quietly = TRUE), "\n")
}, error = function(e) cat("install failed:", conditionMessage(e), "\n"))