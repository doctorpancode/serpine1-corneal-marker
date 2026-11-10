suppressPackageStartupMessages({library(Matrix); library(quadprog); library(scTenifoldKnk)})
D8 <- "D:/A生信分析/虚拟敲除/第8天执行结果"
WT <- as.matrix(readRDS(file.path(D8,"data","D8_WT_network.rds"))); diag(WT) <- 0
genes <- rownames(WT)
Sz <- matrix(rnorm(11*947), 11, 947); colnames(Sz) <- genes
Dmat <- 1e-6*diag(11)
y <- rnorm(947)
cat("dim Sz:", dim(Sz), " dim y:", length(y), "\n")
d <- crossprod(Sz) + Dmat          # PROBLEM: crossprod(Sz) is 947x947, Dmat is 11x11
cat("dim d:", dim(d), " dim Dmat:", dim(Dmat), "\n")
