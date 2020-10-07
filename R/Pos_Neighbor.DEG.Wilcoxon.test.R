argv <- commandArgs(T)
if(length(argv) != 6){stop("Rscript  profile.wilcox.r  [input profile ]  [input name.file ]  [control ] [ case ] [ paired:F/T] [prefix]")}
profile=read.table(argv[1],check.names=F,header=T)
profile<-profile[rowSums(profile)!=0,] 
name=read.table(argv[2],header=T)
config<-table(name[,2])
#sort[,1]=as.character(sort[,1])
#b=subset(profile,select=sort[,1])
data=data.frame()
out=paste(argv[6],".",argv[4],"-vs-",argv[3],".xls",sep="")
#out=paste(argv[4],".xls",sep="")
j=0
all_control=subset(profile,select=as.character(name[name[,2] == argv[3],1]))  
all_case=subset(profile,select=as.character(name[name[,2] == argv[4],1]))
for (i in 1:nrow(profile)){
	#control=as.numeric(profile[i,name[name[,2] == argv[3],1]])
	#case=as.numeric(profile[i,name[name[,2] == argv[4],1]])
	control=as.numeric(all_control[i,])
	case=as.numeric(all_case[i,])
	if(sum(control)+sum(case)>0){
		j=j+1
		data[j,1]=mean(case)
		data[j,2]=mean(control)
		data[j,3]=sd(case)
		data[j,4]=sd(control)
		if(argv[5] == "T"){
			wil=wilcox.test(case,control,paired=T,conf.level=0.95)
			data[j,5]=wil$p.value
			#tmp=t.test(case,control,paired=T,conf.level=0.95)
			#tt=tmp$p.value
		}
		else{
			wil=wilcox.test(case,control,conf.level=0.95)
			data[j,5]=wil$p.value
			#tmp=t.test(case,control,conf.level=0.95)
			#tt=tmp$p.value
		}
		if(data[j,1] == 0 & data[j,2] != 0){
                	data[j,6]="small";
        	}
	        else if(data[j,1] != 0 & data[j,2] == 0){
        	        data[j,6]="big";
	        }
        	else{
                	data[j,6]=data[j,1]/data[j,2];
		}
		
		data[j,7]=data[j,1]-data[j,2]
		rownames(data)[j]=rownames(profile)[i]
       	}
#        else{ exit;
#            }
}
qq<-p.adjust(as.numeric(data[,5]),method="BH")
#qq<-p.adjust(as.numeric(data[,5]),method="holm")
data<-as.data.frame(data)
#rownames(data)<-rownames(profile)
data$qvalue<-qq
colnames(data)=c(paste(argv[4],"Mean",sep="-"),paste(argv[3],"Mean",sep="-"),paste(argv[4],"Sd",sep="-"),paste(argv[3],"Sd",sep="-"),"p-value","multiple",paste(argv[4],"-",argv[3],"(>0:enrich in ",argv[4],"; <0:enrich in ",argv[3],")",sep=""),"qvalue")
write.table(data,file=out,row.names=T,col.names=T,quote=F,sep="\t")
