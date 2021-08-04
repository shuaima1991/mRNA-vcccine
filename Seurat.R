options(stringsAsFactors = F)
rm(list = ls()) 
setwd("")
###################################04.Pre-processing and correction of data###################################
#read data
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
rt=read.table("Inputdata.txt",sep="\t",header=T,check.names=F)


rownames(rt)=rt[,1]
rt=rt[,-1]
rt=round(rt,2)
rt=as.matrix(rt)
rt=log2(rt+1)

exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#The matrix is converted to a Seurat object and the data is filtered
pbmc <- CreateSeuratObject(counts = data,project = "seurat", min.cells = 3, min.features = 50, names.delim = "_",)
#Percentage of mitochondrial genes was calculated using the PercentageFeatureSet function
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
pdf(file="04.featureViolin.pdf",width=10,height=6)          
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
pbmc <- subset(x = pbmc, subset = nFeature_RNA > 50 & percent.mt < 5)    

#Correlation mapping of sequencing depth
pdf(file="04.featureCor.pdf",width=10,height=6)              
plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=1.5)
plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",,pt.size=1.5)
CombinePlots(plots = list(plot1, plot2))
dev.off()

#Standardize the data
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 1500)
#Output characteristic variance graph
top10 <- head(x = VariableFeatures(object = pbmc), 10)
pdf(file="04.featureVar.pdf",width=10,height=6)              #保存基因特征方差图
plot1 <- VariableFeaturePlot(object = pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()





###################################05.PCA###################################
##PCA analysis
pbmc=ScaleData(pbmc)                    
pbmc=RunPCA(object= pbmc,npcs = 20,pc.genes=VariableFeatures(object = pbmc))    

#The genes associated with each PCA component were mapped
pdf(file="05.pcaGene.pdf",width=10,height=8)
VizDimLoadings(object = pbmc, dims = 1:4, reduction = "pca",nfeatures = 20)
dev.off()

#Principal component analysis graph
pdf(file="05.PCA.pdf",width=6.5,height=6)
DimPlot(object = pbmc, reduction = "pca")
dev.off()

#Principal component analysis heat map
pdf(file="05.pcaHeatmap.pdf",width=10,height=8)
DimHeatmap(object = pbmc, dims = 1:4, cells = 500, balanced = TRUE,nfeatures = 30,ncol=2)
dev.off()

#P-value distribution and uniform distribution for each PC
pbmc <- JackStraw(object = pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(object = pbmc, dims = 1:20)
pdf(file="05.pcaJackStraw.pdf",width=8,height=6)
JackStrawPlot(object = pbmc, dims = 1:20)
dev.off()





###################################06.TSNE###################################

pcSelect=20
pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect)                
pbmc <- FindClusters(object = pbmc, resolution = 0.5)                  
pbmc <- RunTSNE(object = pbmc, dims = 1:pcSelect)                      
pdf(file="06.TSNE.pdf",width=6.5,height=6)
TSNEPlot(object = pbmc, do.label = TRUE, pt.size = 2, label = TRUE)   
dev.off()
write.table(pbmc$seurat_clusters,file="06.tsneCluster.txt",quote=F,sep="\t",col.names=F)

logFCfilter=0.5
adjPvalFilter=0.05
pbmc.markers <- FindAllMarkers(object = pbmc,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = logFCfilter)
sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
write.table(sig.markers,file="06.markers.xls",sep="\t",row.names=F,quote=F)
cluster=pbmc@meta.data$orig.ident
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#Heat maps of marker in each cluster are drawn
pdf(file="06.tsneHeatmap.pdf",width=12,height=9)
DoHeatmap(object = pbmc, features = top10$gene) + NoLegend()
dev.off()

#Draw the violin map of marker
pdf(file="06.markerViolin.pdf",width=10,height=6)
VlnPlot(object = pbmc, features = c("IGLL5", "MBOAT1"))
dev.off()

#Draw the scatter diagram of marker in each cluster
pdf(file="06.markerScatter.pdf",width=10,height=6)
FeaturePlot(object = pbmc, features = c("IGLL5", "MBOAT1"),cols = c("green", "red"))
dev.off()

#Draw the bubble diagram of marker in each cluster
pdf(file="06.markerBubble.pdf",width=12,height=6)
cluster10Marker=c("MBOAT1", "NFIB", "TRPS1", "SOX4", "CNN3", "PIM2", "MZB1", "MS4A1", "ELK2AP", "IGLL5")
DotPlot(object = pbmc, features = cluster10Marker)
dev.off()

library(scCATCH)
clu_markers <- findmarkergenes(object = pbmc,
                               species = "Human",
                               cluster = c("IS1", "IS2", "IS3", "IS4"),
                               match_CellMatch = TRUE,
                               cancer = NULL,
                               tissue = c('Brain'))

clu_markers <- findmarkergenes(object = pbmc,
                               species = "Human",
                               cluster = 'All',
                               match_CellMatch = FALSE,
                               cancer = NULL,
                               tissue = NULL,
                               cell_min_pct = 0.25,
                               logfc = 0.25,
                               pvalue = 0.05)

clu_ann <- scCATCH(object = clu_markers$clu_markers,
                   species = "Human",
                   cancer = NULL,
                   tissue = 'Brain')


write.table(clu_ann,"clu_ann.txt",sep = "\t")

write.table(clu_markers,"clu_markers.txt",sep = "\t")




###################################07.Annotated cell type###################################




library(SingleR)
counts<-pbmc@assays$RNA@counts
clusters<-pbmc@meta.data$seurat_clusters
ann=pbmc@meta.data$orig.ident
singler = CreateSinglerObject(counts, annot = ann, "pbmc", min.genes = 0,
  species = "Human", citation = "",
  ref.list = list(), normalize.gene.length = F, variable.genes = "de",
  fine.tune = F, do.signatures = T, clusters = clusters, do.main.types = T,
  reduce.file.size = T, numCores = 1)
singler$seurat = pbmc
singler$meta.data$xy = pbmc@reductions$tsne@cell.embeddings
clusterAnn=singler$singler[[2]]$SingleR.clusters.main$labels
write.table(clusterAnn,file="07.clusterAnn.txt",quote=F,sep="\t",col.names=F)
write.table(singler$other,file="07.cellAnn.txt",quote=F,sep="\t",col.names=F)

#Prepare documents needed for Monocle analysis
monocle.matrix=as.matrix(pbmc@assays$RNA@data)
monocle.matrix=cbind(id=row.names(monocle.matrix),monocle.matrix)
write.table(monocle.matrix,file="07.monocleMatrix.txt",quote=F,sep="\t",row.names=F)
monocle.sample=as.matrix(pbmc@meta.data)
monocle.sample=cbind(id=row.names(monocle.sample),monocle.sample)
write.table(monocle.sample,file="07.monocleSample.txt",quote=F,sep="\t",row.names=F)
monocle.geneAnn=data.frame(gene_short_name = row.names(monocle.matrix), row.names = row.names(monocle.matrix))
monocle.geneAnn=cbind(id=row.names(monocle.geneAnn),monocle.geneAnn)
write.table(monocle.geneAnn,file="07.monocleGene.txt",quote=F,sep="\t",row.names=F)
write.table(singler$other,file="07.monocleClusterAnn.txt",quote=F,sep="\t",col.names=F)
write.table(sig.markers,file="07.monocleMarkers.txt",sep="\t",row.names=F,quote=F)