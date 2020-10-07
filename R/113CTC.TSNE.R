library(ggplot2)
library(Rtsne)

set.seed(54)

exp <- read.table("113CTC_tumor_cellline.m1n13.exp.matrix", header = T)
anno <- read.table("anno.txt", header = T)

exp <- t(log2(exp + 1))


tsne <- Rtsne(test, check_duplicates = FALSE,
              pca = TRUE, perplexity=8, theta=0, dims=2)

embedding <- as.data.frame(tsne$Y)
embedding$Group <- as.factor(anno$Group[anno$Group=="CTC"])
embedding$PID <- as.factor(embedding$PID, levels=paste0("P", 1:10))
write.table(embedding, "CTC.tsne.anno.txt", quote=F, sep="\t", row.names = F, col.names = T)

pdf("113CTC.TSNE.pdf", 5.5, 4.5)
ggplot(embedding, aes(x=V1, y=V2, color=PID),alpha=0.30) +
  geom_point(size=1.5) + 
  #scale_colour_brewer(palette="Spectral") +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "red", "black")) +
  theme_bw() + 
  theme(panel.grid =element_blank()) +
  labs(x= "tSNE.1",y = "tSNE.2")
dev.off()