rm(list = ls())
library(ggplot2)
library(ggpubr)
library(ggExtra)
setwd("")            
pFilter=0.01       

#read data
exp=read.table("singleGeneExp.txt", header=T,sep="\t",check.names=F)
exp=exp[!duplicated(exp$id), ]
rownames(exp)=exp[,1]
exp=exp[,-1]
gene=colnames(exp)[1]

TME=read.table("estimateScores1.txt", header=T,sep="\t",check.names=F)
TME=TME[!duplicated(TME$ID), ]
rownames(TME)=TME[,1]
TME=TME[,-1]

#Sample intersection
sameSample=intersect(row.names(TME),row.names(exp))
TME=TME[sameSample,]
exp=exp[sameSample,]

# correlation test 
outTab=data.frame()
#circulation 
for(i in as.factor(exp[,"CancerType"])){
    exp1=exp[(exp[,"CancerType"]==i),]
    TME1=TME[(TME[,"CancerType"]==i),]
    y=as.numeric(exp1[,1])
    outVector=data.frame(i,gene)
	
	for(j in colnames(TME1)[1:10]){
		x=as.numeric(TME1[,j])
		df1=as.data.frame(cbind(x,y))
		corT=cor.test(x,y,method="spearman")
		cor=corT$estimate
		pValue=corT$p.value
		outVector=cbind(outVector,pValue)
		p1=ggplot(df1, aes(x, y)) + 
			xlab(j)+ylab(gene)+
			ggtitle(paste0("Cancer: ",i))+theme(title=element_text(size=10))+
		    geom_point()+ geom_smooth(method="lm") + theme_bw()+
		    stat_cor(method = 'spearman', aes(x =x, y =y))
	    p2=ggMarginal(p1, type = "density", xparams = list(fill = "orange"),yparams = list(fill = "blue"))
		if(pValue<pFilter){
			pdf(file=paste0("estimateCor.",i,"_",j,".pdf"),width=5,height=5)
			print(p2)
			dev.off()
		}
	}
	outTab=rbind(outTab,outVector)
}
colNames=c("CancerType","Gene",colnames(TME)[1:10])
colnames(outTab)=colNames
write.table(outTab,file="estimateCor.result.txt",sep="\t",row.names=F,quote=F)
