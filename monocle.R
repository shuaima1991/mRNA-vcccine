options(stringsAsFactors = F)
rm(list = ls()) 
setwd("")
library(monocle)
monocle.matrix=read.table("07.monocleMatrix.txt",sep="\t",header=T,row.names=1,check.names=F)
monocle.sample=read.table("07.monocleSample.txt",sep="\t",header=T,row.names=1,check.names=F)
monocle.geneAnn=read.table("07.monocleGene.txt",sep="\t",header=T,row.names=1,check.names=F)
marker=read.table("07.monocleMarkers.txt",sep="\t",header=T,check.names=F)

#Convert the Seurat results into the cell matrix, cell annotation table and gene annotation table required by Monocle
data <- as(as.matrix(monocle.matrix), 'sparseMatrix')
pd<-new("AnnotatedDataFrame", data = monocle.sample)
fd<-new("AnnotatedDataFrame", data = monocle.geneAnn)
cds <- newCellDataSet(data, phenoData = pd, featureData = fd)

#Rename one of the columns
names(pData(cds))[names(pData(cds))=="seurat_clusters"]="Cluster"
pData(cds)[,"Cluster"]=paste0("cluster",pData(cds)[,"Cluster"])

#Add cell clustering data
clusterRt=read.table("clu_ann1.txt",header=F,sep="\t",check.names=F)
clusterAnn=as.character(clusterRt[,2])
names(clusterAnn)=paste0("cluster",clusterRt[,1])

#After the comment is completed, just put PDATA (CD)$to select your own category
pData(cds)$cell_type2 <- plyr::revalue(as.character(pData(cds)$Cluster),clusterAnn)

#Pseudo-time analysis process
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds <- setOrderingFilter(cds, marker$gene)
plot_ordering_genes(cds)
cds <- reduceDimension(cds, max_components = 2,reduction_method = 'DDRTree')
cds <- orderCells(cds)
pdf(file="cluster.trajectory1.pdf",width=6.5,height=6)
plot_cell_trajectory(cds, color_by = "Cluster")
dev.off()
pdf(file="cellType.trajectory1.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "cell_type2")
dev.off()
