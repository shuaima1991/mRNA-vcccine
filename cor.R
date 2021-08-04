rm(list = ls())
gene="riskScore"           
geneRT=read.table("gene.txt",sep="\t",header=F)          
files='stemness4-LGG.txt'
data=read.table(files, header=T,sep="\t",check.names=F,row.names = 1)
CancerType=levels(data$CancerType)
sameGenes=sameGenes=intersect(as.vector(geneRT[,1]),colnames(data))
outTab=matrix(1,31,4)

corTab=matrix(1,31,4)

for(i in 1:31){
  data1=data[data$CancerType==CancerType[i],]
  x=as.numeric(data1[,gene])
  
  for(j in 1:4){
    y=as.numeric(data1[,j])
    corT=cor.test(x,y)
    cor=corT$estimate
    pValue=corT$p.value
    
    outTab[i,j]=pValue
    corTab[i,j]=cor
  }}

colnames(outTab)=sameGenes
row.names(outTab)=CancerType
write.table(outTab,file="geneCor.pvalue.txt",sep="\t",row.names=T,quote=F)
colnames(corTab)=sameGenes
row.names(corTab)=CancerType
write.table(corTab,file="geneCor.cor.txt",sep="\t",row.names=T,quote=F)

