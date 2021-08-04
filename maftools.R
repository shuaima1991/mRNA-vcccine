setwd("")
library(maftools)
luad <- read.maf(maf="TCGA.GBM.LGG.maf")

# Extracting gender-specific "Tumor_Sample_Barcode" from clinical data
clin <- read.table("TCGA-GBM.LGG_phenotype .tsv", header=T, sep="\t")
clin.IS1 <- subset(clin, Type=="IS1")$Tumor_Sample_Barcode
clin.IS2 <- subset(clin, Type=="IS2")$Tumor_Sample_Barcode
clin.IS3 <- subset(clin, Type=="IS3")$Tumor_Sample_Barcode
clin.IS4 <- subset(clin, Type=="IS4")$Tumor_Sample_Barcode
# Use SubsetMAF to build male and female MAF objects
luad.IS1 <- subsetMaf(maf=luad, tsb=clin.IS1, isTCGA=TRUE)
luad.IS2 <- subsetMaf(maf=luad, tsb=clin.IS2, isTCGA=TRUE)
luad.IS3 <- subsetMaf(maf=luad, tsb=clin.IS3, isTCGA=TRUE)
luad.IS4 <- subsetMaf(maf=luad, tsb=clin.IS4, isTCGA=TRUE)

#drive genes
shnsc.sigIS1 = oncodrive(maf = luad.IS1, AACol = 'HGVSp_Short', minMut = 5, pvalMethod = 'zscore')
plotOncodrive(res = shnsc.sigIS1, fdrCutOff = 0.01, useFraction = TRUE,labelSize = 1)
shnsc.sigIS2 = oncodrive(maf = luad.IS2, AACol = 'HGVSp_Short', minMut = 5, pvalMethod = 'zscore')
plotOncodrive(res = shnsc.sigIS2, fdrCutOff = 0.01, useFraction = TRUE,labelSize = 1)
shnsc.sigIS3 = oncodrive(maf = luad.IS3, AACol = 'HGVSp_Short', minMut = 5, pvalMethod = 'zscore')
plotOncodrive(res = shnsc.sigIS3, fdrCutOff = 0.01, useFraction = TRUE,labelSize = 1)
shnsc.sigIS4 = oncodrive(maf = luad.IS4, AACol = 'HGVSp_Short', minMut = 5, pvalMethod = 'zscore')
plotOncodrive(res = shnsc.sigIS4, fdrCutOff = 0.01, useFraction = TRUE,labelSize = 1)


#Use the color of RColorBrewer here, but you can use any color
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)
#View the color corresponding to the variation type
print(vc_cols)
#>   Frame_Shift_Del Missense_Mutation Nonsense_Mutation         Multi_Hit 
#>         "#A6CEE3"         "#1F78B4"         "#B2DF8A"         "#33A02C" 
#>   Frame_Shift_Ins      In_Frame_Ins       Splice_Site      In_Frame_Del 
#>         "#FB9A99"         "#E31A1C"         "#FDBF6F"         "#FF7F00"

oncoplot(luad.IS1, colors = vc_cols, top = 20)
oncoplot(luad.IS2, colors = vc_cols, top = 20)
oncoplot(luad.IS3, colors = vc_cols, top = 20)
oncoplot(luad.IS4, colors = vc_cols, top = 20)
output <- somaticInteractions(maf=luad.IS1, top=25, pvalue=c(0.05, 0.01))
output <- somaticInteractions(maf=luad.IS2, top=25, pvalue=c(0.05, 0.01))
output <- somaticInteractions(maf=luad.IS3, top=25, pvalue=c(0.05, 0.01))
output <- somaticInteractions(maf=luad.IS4, top=25, pvalue=c(0.05, 0.01))

#VAF
plotVaf(maf = luad.IS1,width = 100,height = 10)
plotVaf(maf = luad.IS2,width = 100,height = 10)
plotVaf(maf = luad.IS3,width = 100,height = 10)
plotVaf(maf = luad.IS4,width = 100,height = 10)
luad.vaf <- vafCompare(m1 = luad.IS2,m2 = luad.IS3)


# 
fvsm <- mafCompare(m1=luad.IS1, m2=luad.IS2, m3=luad.IS3, m4=luad.IS4, m1Name="IS1", m2Name="IS2",m3Name="IS3",m4Name="IS4", minMut=5)
# 结果保存到文件"female_vs_male.tsv"
write.table(fvsm$results, file="female_vs_male.tsv", quote=FALSE, row.names=FALSE, sep="\t")


#forest map 
fvsm1 <- mafCompare(m1=luad.IS1, m2=luad.IS2, m1Name="IS1", m2Name="IS2", minMut=5)
forestPlot(mafCompareRes=fvsm1, pVal=0.001, color=c("maroon", "royalblue"), geneFontSize=0.8)
fvsm2 <- mafCompare(m1=luad.IS2, m2=luad.IS3, m1Name="IS2", m2Name="IS3", minMut=5)
forestPlot(mafCompareRes=fvsm2, pVal=0.001, color=c("maroon", "royalblue"), geneFontSize=0.8)
fvsm3 <- mafCompare(m1=luad.IS3, m2=luad.IS4, m1Name="IS3", m2Name="IS4", minMut=5)
forestPlot(mafCompareRes=fvsm3, pVal=0.001, color=c("maroon", "royalblue"), geneFontSize=0.8)
fvsm4 <- mafCompare(m1=luad.IS1, m2=luad.IS3, m1Name="IS1", m2Name="IS3", minMut=5)
forestPlot(mafCompareRes=fvsm3, pVal=0.001, color=c("maroon", "royalblue"), geneFontSize=0.8)
