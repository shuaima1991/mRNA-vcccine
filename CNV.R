
setwd("D:/Rcode/文章思路/细胞焦亡与抗肿瘤免疫/4分组突变/CNV/Low immunityeraly")

#读取gistic文件
library(maftools)
luad.gistic <- readGistic(gisticAllLesionsFile="all_lesions.conf_90.txt.txt", gisticAmpGenesFile="amp_genes.conf_90.txt.txt", gisticDelGenesFile="del_genes.conf_90.txt.txt", gisticScoresFile="scores (3).gistic", isTCGA=TRUE)

#染色体图
gisticChromPlot(gistic=luad.gistic, markBands="all")


#气泡图
gisticBubblePlot(gistic=luad.gistic)


#瀑布图
luad <- read.maf(maf="TCGA.LUAD.maf", clinicalData="LUAD-TP.samplefeatures.txt")
gisticOncoPlot(gistic=luad.gistic, clinicalData=getClinicalData(x=luad), clinicalFeatures="CLI_gender", sortByAnnotation=TRUE, top=10)