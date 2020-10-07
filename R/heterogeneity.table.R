# Usage: Rscript heterogeneity_pos_plot.R <XX.exp>
args = commandArgs(T)
library(bootstrap)
exp <- read.table(args[1],head=T,check.names=F)
###group
group = c("HV","PA","PV","PoV","HV|PA|PV|PoV")
for(i in 1:length(group)){
  test = exp[,grep(group[i],colnames(exp))]
  if(length(test)==0||length(ncol(test))==0){
    M=s=0
    c=0
 }
  else{
#########filter genes
  isexpr <- rowSums(test>1)>ncol(test)*0.5 ##for a gene: the number of samples with fpkm>1 must exceed the half.
  test <- test[isexpr,]
  LOG = log2(test+1)
  LOG = transform(LOG,var=apply(LOG,1,var))
  LOG = LOG[order(LOG$var,decreasing=T),]
  len = length(colnames(LOG))
  LOG = LOG[,-len]
  top2000 = head(LOG,n=2000)
##########M calculation
  c = round(cor(top2000[,1:(len-1)],top2000[,1:(len-1)]),4)
  tmp =c # tmp = atanh(c)
  M=mean(tmp[lower.tri(tmp)])
##########sd calculation
  theta <- function(x){mean(x)}
  results <- jackknife(tmp[lower.tri(tmp)],theta)
  if (length(grep("TRUE",is.na(results$jack.values)))>0){
  #if(is.na(results$jack.values)){
	s=0
  }else{
  s=sd(results$jack.values)*qnorm(0.975)
  }  
}
  if(i==1) {
  m = tanh(M)
  s1 = tanh(M)-s
  s2 = tanh(M)+s
  }
  else {
  m = c(m,tanh(M))
  s1 = c(s1,tanh(M)-s)
  s2 = c(s2,tanh(M)+s)
  }
cor = c[lower.tri(c)]
m_sd = data.frame(m = m, s1 = s1, s2 =s2)
write.table(cor,file=paste(args[2],group[i],"cor.table.txt",sep="."),row.names=F,col.names=F,sep="\t",quote=F)
write.table(m_sd,file=paste(args[2],group[i],"m.table.txt",sep="."),row.names=F,col.names=F,sep="\t",quote=F)
}

