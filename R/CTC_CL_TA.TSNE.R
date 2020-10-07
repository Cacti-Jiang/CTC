library(ggplot2)
library(Rtsne)

set.seed(13)

exp <- read.table("113CTC_tumor_cellline.m1n13.exp.matrix", header = T)
anno <- read.table("anno.txt", header = T)
exp <- t(log2(exp + 1))

tsne <- Rtsne(test, check_duplicates = FALSE,
              pca = TRUE, perplexity=30, theta=0, dims=2)

embedding <- as.data.frame(tsne$Y)
embedding$Group <- as.factor(anno$Group)

pdf("CTC_CL_TA.TSNE.pdf", 5.5, 4.5)
ggplot(embedding, aes(x=V1, y=V2, color=Group),alpha=0.30) +
  geom_point(size=1.5) + 
  scale_color_manual(values=c(CTC="blue", "cell_line"="green", "primary_tumor"="red")) +
  theme_bw() + 
  theme(panel.grid =element_blank()) +
  labs(x= "tSNE.1",y = "tSNE.2")
dev.off()