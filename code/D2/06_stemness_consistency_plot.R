suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
out <- "D:/A生信分析/虚拟敲除/第2天执行结果"
d2 <- fread(file.path(out,"data","stemness_score.tsv"))
d2[, binary := factor(binary, levels=c("Normal","Lesion"))]
p1 <- ggplot(d2, aes(binary, StemnessScore, fill=binary)) +
  geom_boxplot(width=.55, outlier.shape=NA, alpha=.7) +
  geom_jitter(width=.1,size=2) +
  labs(x=NULL,y="StemnessScore (ssGSEA)",subtitle="GSE241715 (n=21)  P=0.0056") +
  scale_fill_manual(values=c(Normal="#4C72B0",Lesion="#C44E52")) + theme_bw(12)+theme(legend.position="none")

d3 <- fread(file.path(out,"data","GSE204791_stemness_score.tsv"))
d3[, grp := factor(paste(compartment,group,sep="_"), levels=c("Epithelium_Control","Epithelium_KC","Stroma_Control","Stroma_KC"))]
p2 <- ggplot(d3, aes(grp,StemnessScore,fill=group)) +
  geom_boxplot(width=.55,outlier.shape=NA,alpha=.7) +
  geom_jitter(width=.1,size=2) +
  labs(x=NULL,y="StemnessScore (ssGSEA)",subtitle="GSE204791 (n=31)  both P>0.05, same direction") +
  scale_fill_manual(values=c(Control="#4C72B0",KC="#C44E52")) + theme_bw(12) +
  theme(axis.text.x=element_text(angle=25,hjust=1),legend.position="none")
p <- p1 | p2
ggsave(file.path(out,"figures","D2_Fig9_stemness_consistency.pdf"), p, width=10, height=4.4)
ggsave(file.path(out,"figures","D2_Fig9_stemness_consistency.png"), p, width=10, height=4.4, dpi=300)
cat("ok\n")