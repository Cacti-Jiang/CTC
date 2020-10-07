argv <- commandArgs(T)
if(length(argv) != 3){stop("Rscript  kruskal.test.r  [input gene.exp.matrix ] [input group.info ]  [output test.txt ] ")}
data <- read.table(argv[1],check.names=F,head=T)
data <- data[rowSums(data)!=0,]
config <- read.table(argv[2],check.names=F,head=F)
table(config[,2])
configg<-table(config[,2])
A<-which(config[,2]==rownames(configg)[1])
B<-which(config[,2]==rownames(configg)[2])
C<-which(config[,2]==rownames(configg)[3])
D<-which(config[,2]==rownames(configg)[4])
res<-matrix(0,dim(data)[1],5)
for(i in 1:dim(data)[1]){
        gene<-data[i,]
		kk<-as.data.frame(t(gene))
		kk$info<-config[,2]
        res[i,1]<-kruskal.test(kk[,1],kk[,2])$p.value
		rank_all<-as.numeric(rank(gene))
        res[i,2]<-mean(rank_all[A])
		res[i,3]<-mean(rank_all[B])
		res[i,4]<-mean(rank_all[C])
		res[i,5]<-mean(rank_all[D])
}
#qq<-p.adjust(as.numeric(res[,1]),method="BH")
qq<-p.adjust(as.numeric(res[,1]),method="holm")
ress<-as.data.frame(res)
rownames(ress)<-rownames(data)
colnames(ress)<-c("p.value",rownames(configg))
ress$qvalue<-qq
write.table(ress,file=argv[3],sep="\t",quote=F,col.names=NA)

