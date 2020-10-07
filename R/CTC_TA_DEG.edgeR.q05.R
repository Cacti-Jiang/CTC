#Usage: Rscript 2.diff.R sample.lst count.matrix.txt  [prefix]
args <- commandArgs(T)
library(edgeR)

sample<-read.table(args[1],header=T)
data <- read.table(args[2], check.names=FALSE, stringsAsFactors=FALSE)
prefix <-args[3]
data <- data[, rownames(sample)]

group<-factor(sample$Group)
dge <- DGEList(counts=data,group=group, genes=rownames(data))
keep <- rowSums(cpm(dge)>1) >= 2
DE <- dge[keep,]
######### normalize
DE$samples$lib.size <- colSums(DE$counts)
DE <- calcNormFactors(DE)
design <- model.matrix(~group)
rownames(design) <- colnames(DE)
######### Estimating the dispersion
DE <- estimateGLMCommonDisp(DE, design, verbose=TRUE)
DE <- estimateGLMTrendedDisp(DE, design)
DE <- estimateGLMTagwiseDisp(DE, design)
######### Differential expression
fit <- glmFit(DE, design)  ##### Fit genewise glms
lrt <- glmLRT(fit)         ##### likelihood ratio tests for tumour vs normal tissue differences
summary(de<-decideTestsDGE(lrt, p=0.05, adjust="BH")) ######## adjust P value
edge <- as.data.frame(topTags(lrt, n=50000))
#count <- cbind(edge,DE$counts[rownames(edge),])

########## filter FDR <= 0.05  &&   up-regulate/down-regulate
edge2fold_up <- edge[edge$logFC >=1 & edge$FDR <= 0.05,]
#count_up<-cbind(edge2fold_up,DE$counts[rownames(edge2fold_up),])
write.table(edge2fold_up,file=paste(prefix,"up.q05.txt",sep="."),col.names=T,row.names=F,sep="\t",quote=F)

edge2fold_down <- edge[edge$logFC <= -1 & edge$FDR <= 0.05,]
#count_down<-cbind(edge2fold_down,DE$counts[rownames(edge2fold_down),])
write.table(edge2fold_down,file=paste(prefix,"down.q05.txt",sep="."),col.names=T,row.names=F,sep="\t",quote=F)

write.table(edge, "allgene.result.txt", col.names=T,row.names=F,sep="\t",quote=F)
