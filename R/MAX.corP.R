library(plyr)
exp = read.table("CTC.exp.txt", header = T, check.names = F, stringsAsFactors = T)
anno <- read.table("CTC.pos.txt", header = T)

gene <- c("MAX", "CCL5", "MYC", "MAPK14", "MAPK11", "MAPK12", "MAPK13")
mat <- exp[rownames(exp)%in%gene, ]
mat <- log2(mat + 1)

max <- as.numeric(mat["MAX", ])
M_MAX <- M_MAX.p <- rep(NA, 6)  
for (i in 2:nrow(mat)){
  M_MAX[i] <- cor(as.numeric(mat[i, ]), max)
  M_MAX.p[i] <- cor.test(as.numeric(mat[i, ]), max)$p.value
}
cormat <- cbind(M_MAX, M_MAX.p)
rownames(cormat) <- rownames(mat)
write.table(cormat, file = "MAX.P.txt", row.names = T, col.names = NA, sep = "\t", quote = F)

##plot
pdf("MAX.corP.exp.pdf")
par(mfrow=c(3,3), mgp = c(1.5, 0.5, 0), mar = c(3, 3, 3, 3), oma = c(4, 2, 3, 0.5))
for (i in 2:7){
  plot(max, mat[rownames(cormat)[i], ], xlab = "MAX", ylab = rownames(cormat)[i], 
       col = rgb(1, 0, 0, alpha = 0.6), pch = 20)
  mtext(side=3,text=paste("r = ", round(cormat[i, ][1], 3), "P = ", round(cormat[i, ][2], 3), "N = 113", sep = " ")
        ,line=1.2, cex = 0.6)
  mtext(side=3,text="method : Pearson", line=0.5, cex = 0.6)
  abline(lm(as.numeric(mat[rownames(cormat)[i], ])~max))
}
dev.off()
