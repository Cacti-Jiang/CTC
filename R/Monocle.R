## Usage: Rscript Monocle.R <FPKM_matrix.txt> <Cell_sample_sheet.txt> <gene_annotations.txt> <out.prefix> <num_path> 

library("monocle")
library("reshape2")

args <- commandArgs(T)

## load data
expr_matrix <- read.table(args[1],header=T,check.names=F)
sample_sheet <- read.table(args[2],header=T)
sample_sheet <- as.data.frame(t(sample_sheet))
gene_annotation <- read.table(args[3],header=T)
"load is done"


## creat dataset and retain genes expressed in at least 10% cells with FPKM  >1
pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)
SAMP <- newCellDataSet(as.matrix(expr_matrix), phenoData = pd, featureData = fd)
SAMP <- detectGenes(SAMP, min_expr = 1)
cell_len<-dim(expr_matrix)[2]
expressed_genes <- row.names(subset(fData(SAMP), num_cells_expressed > cell_len*0.1))
#write.table(SAMP[expressed_genes,], sep="\t", file=paste(args[4],"diff.test2.txt",sep="."),quote=FALSE)
"dataset is created"


## find differentially expressed genes
diff_test_res <- differentialGeneTest(SAMP[expressed_genes,],fullModelFormulaStr="expression~Group", cores=3)
write.table(diff_test_res, sep="\t", file=paste(args[4],"diff.test.txt",sep="."),quote=FALSE)
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
ordering_genes <- intersect(ordering_genes, expressed_genes)
"differentially expressed genes is detected"


## order cells by progress
SAMP <- setOrderingFilter(SAMP, ordering_genes)
SAMP <- reduceDimension(SAMP, use_irlba=FALSE) ## Reducing to independent components so Monocle will work better
save(SAMP, file=paste(args[4],"ori.Rdata",sep="."))
SAMP <- orderCells(SAMP, num_paths=as.numeric(args[5]), reverse=TRUE) ##  num_paths allows Monocle to assign cells to one of several alternative fates
pdf(file=paste(args[4],"ordering.pdf",sep="."))
#plot_spanning_tree(SAMP,show_cell_names =F,cell_name_size = 1.8,color_by="State")
plot_spanning_tree(SAMP,show_backbone=F,show_cell_names =F,cell_name_size = 1.8,color_by="Group")
dev.off()
"ordering is done"

###print order of pseudotime express
vara<-pData(SAMP)$Pseudotime
names(vara)<-rownames(pData(SAMP))
pseu_exprs<-t(expr_matrix)
pseu_exprs<-pseu_exprs[names(sort(vara)),ordering_genes]
pseu_exprs<-t(pseu_exprs)
write.table(pseu_exprs,file=paste(args[4],"pseudotime.exprs.txt",sep="."),sep="\t",quote=F)