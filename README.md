# CRC_ovarian-metastasis

rm(list = ls())
gc()
options(stringsAsFactors = F)
library(dplyr)
library(stringr)
library(maftools)
setwd('0330all54-maf')
var_maf = read.maf(maf ="0330all54.maf")
#maf2 = read.maf(maf ="0330all63.maf")
sam="0501All54"
maf="0330all54.maf"
laml.maf = read.maf(maf = maf, verbose = FALSE)
#laml.maf@data <- gsub('shl_660va','shl_66ova',laml.maf@data)
laml.maf@data[laml.maf@data$Tumor_Sample_Barcode=='shl_660va',1]<-'shl_66ova'

samp=getSampleSummary(laml.maf)
samp
write.table(samp, file="0501all54_summary2.txt", quote= F,sep="\t", row.names=F)
laml.mutsig=read.table("Dndcv_musig_q0.05.txt",header=T,dec=".",fill=T)
head(laml.mutsig)
write.table(laml.mutsig,file=paste0(sam,"DNDcv-sig.txt"),sep="\t",row.names=F)

laml.clin2 = read.table("0330all54_clinical.txt",header=T,dec=".",fill=T)
samp_clin2<-unique(laml_clin2$Tumor_Sample_Barcode)
samp_clin2
write.table(laml_clin2, file="0502all54_clinical.txt", quote= F, sep="\t", row.names=F)

laml_maf3 = subsetMaf(maf = laml.maf, query = "Tumor_Sample_Barcode != 'xyx_92per'")
laml_maf3@data[laml_maf3@data$Tumor_Sample_Barcode=='xyy_87nd_',1]<-'xyy_87nod'
samp2<-unique(laml_maf3@data$Tumor_Sample_Barcode)
samp2
save(laml_maf3,file='0502All54snv.Rdata')
laml54<-get(load('0502All54snv.Rdata'))
head(laml54@data)
write.mafSummary(laml_maf3, basename = "202205all54", compress = FALSE)
maf2='202205all54_maftools.maf'

#CNV
setwd('20220502_all54_maftools')
sam5="0502all54"
all.lesions="all_lesions.conf_75.txt"
amp.genes = "amp_genes.conf_75.txt"
del.genes = "del_genes.conf_75.txt"
scores.gis2 = "scores.gistic"

laml.plus.gistic2 = read.maf(
  maf=maf2,
  gisticAllLesionsFile = all.lesions,
  gisticAmpGenesFile = amp.genes,
  gisticDelGenesFile = del.genes,
  gisticScoresFile = scores.gis2,
  isTCGA = FALSE,
  verbose = FALSE,
  clinicalData = laml_clin2
)
head(laml.plus.gistic2@data$Tumor_Sample_Barcode)
laml.plus.gistic2 = subsetMaf(maf = laml.plus.gistic2, query = "Tumor_Sample_Barcode != 'xyy_87nd_'")
head(laml.plus.gistic2@data$Tumor_Sample_Barcode)

save(laml.plus.gistic2,file='0502all54-cnv-snv-clinical.Rdata')
laml.gistic2 = readGistic(gisticAllLesionsFile = all.lesions, gisticAmpGenesFile = amp.genes, gisticDelGenesFile = del.genes, gisticScoresFile = scores.gis2, isTCGA = F)
pdf(paste0(sam5, "_gistic1.pdf"),width=8,height=8)
gisticChromPlot(gistic = laml.gistic2, markBands = "all")
dev.off()
pdf(paste0(sam5, "_gistic_bubble.pdf"),width=8,height=8)
gisticBubblePlot(gistic = laml.gistic2)
dev.off()
pdf(paste0(sam5, "_gistic_onco10.pdf"),width=28,height=8)
gisticOncoPlot(gistic = laml.gistic2, clinicalData = getClinicalData(x = laml.plus.gistic2), clinicalFeatures = 'metaSite', sortByAnnotation = TRUE, top = 10, fontSize=1.2, gene_mar=10, legendFontSize=1.6, showTumorSampleBarcodes=T, SampleNamefontSize=1)
dev.off()

pdf(paste0(sam5, "_TMB.pdf"),width=8,height=8)
tmb(maf = laml.plus.gistic2, captureSize = 10)
dev.off()
   
##Compare mutation load against TCGA cohorts
tcgaCompare uses mutation load from TCGA MC3 for comparing muttaion burden against 33 TCGA cohorts. Plot generated is similar to the one described in Alexandrov et al 5.
samp2="All54"
pdf(paste0(sam5,"_TMB-COAD", ".pdf"),width=8,height=8)
laml.mutload = tcgaCompare(maf = laml.plus.gistic2, cohortName = samp2, logscale = TRUE, capture_size = 10)
dev.off()
write.table(laml.mutload$median_mutation_burden, file=paste0(samp2, "_TMB_COAD_diff.csv"))
save(laml.mutload, file=paste0(samp2, "_TMB_coad_diff.Rdata"))


#Shows sample summry.
getSampleSummary(laml.plus.gistic2)
#Shows gene summary.
getGeneSummary(laml.plus.gistic2)
#shows clinical data associated with samples
getClinicalData(laml.plus.gistic2)
#Shows all fields in MAF
getFields(laml.plus.gistic2)
#Writes maf summary to an output file with basename laml.
#write.mafSummary(maf = laml.plus.gistic, basename = 'laml.plus.gistic_summary')
pdf(paste0(sam5, "_summary.pdf"),width=8,height=8)
plotmafSummary(maf = laml.plus.gistic2, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE, showBarcodes=TRUE, top=20)
dev.off()


sample<-unique(laml.plus.gistic2@data$Tumor_Sample_Barcode)
pdf(paste0(sam5,"_maftool_everything50.pdf"), width=18, height=18)
oncoplot(
  maf = laml.plus.gistic2,
  top=50,  #draw_titv = TRUE,
  pathways = "auto",
  clinicalFeatures = c('Age35','Age50', 'Path', 'Location', 'Lymph', 'MSS', 'OS12', 'OS24'),
  #sortByAnnotation = TRUE,
  additionalFeature = c("Tumor_Seq_Allele2", "C"),
  #leftBarData = aml_genes_vaf,
  #leftBarLims = c(0, 100),
  rightBarData = laml.mutsig,
  #rightBarData = laml.mutsig[,.(gene, q)],
titleFontSize=1.4,
legendFontSize=1.4,
removeNonMutated=FALSE,
fontSize = 0.7,
showTumorSampleBarcodes=T,
SampleNamefontSize=0.9,
sampleOrder=sample,
annotationFontSize=1.2,
writeMatrix = T
)
dev.off()

##all54_21-all54 vs Coad

driver_genes23= c("KRAS", "NRAS", "TP53", "APC", "BRAF", "RNF43", "PCDHB12", "ACVR2A", "ZNF160", "ZNF716", "STOML1", "SMIM3", "NLGN1", "DMD", "LRP2", "FAT4", "ARID1A", "NCOR1", "DVL2", "TLE3", "MNT", "MUC16", "PRKDC")
pdf(paste0("0709all54_Coad","_drivergenes23_oncoplot.pdf"), width=18, height=18)
 oncoplot(
  maf = laml.plus.gistic2,
   #top=200,
   genes=driver_genes23,
   draw_titv = F,
   pathways = NULL,
   gene_mar=15,
   clinicalFeatures = c('Age35','Age50', 'Path', 'Location', 'Lymph', 'MSS', 'OS24'),
   #sortByAnnotation = TRUE,
  additionalFeature = c("Tumor_Seq_Allele2", "C"),
   #leftBarData = aml_genes_vaf,
   #leftBarLims = c(0, 100),
   #rightBarData = laml.mutsig,
   #rightBarData = laml.mutsig[,.(gene, q)],
 titleFontSize=1.4,
 legendFontSize=1.4,
 removeNonMutated=FALSE,
 fontSize = 0.9,
 showTumorSampleBarcodes=T,
 SampleNamefontSize=0.9,
 sampleOrder=sample,
 annotationFontSize=1.2,
 writeMatrix = T
 )
 dev.off()

pdf(paste0("0709coad","_drivergenes23_oncoplot.pdf"), width=18, height=18)
 oncoplot(
 maf = coad.maf4,
   #top=200,
    genes=driver_genes23,
   draw_titv = F,
    pathways = NULL,
    gene_mar=12,
    #clinicalFeatures = c('Age35','Age50', 'Path', 'Location', 'Lymph', 'MSS', 'OS24'),
    #sortByAnnotation = TRUE,
  additionalFeature = c("Tumor_Seq_Allele2", "C"),
    #leftBarData = aml_genes_vaf,
   #leftBarLims = c(0, 100),
    #rightBarData = laml.mutsig,
   #rightBarData = laml.mutsig[,.(gene, q)],
  titleFontSize=1.4,
  legendFontSize=1.4,
  #removeNonMutated=FALSE,
  fontSize = 0.6,
  #showTumorSampleBarcodes=F,
  #SampleNamefontSize=0.9,
  #sampleOrder=sample,
  annotationFontSize=1.2,
  writeMatrix = F
  )
  dev.off()


#all54-coad-differ
os01 <- mafCompare(m1 = laml.plus.gistic3, m2 = coad.maf2, m1Name = 'All54', m2Name = 'COAD460', minMut = 20, useCNV=TRUE)
print(os01)
write.table(os01$result,file=paste0("all54-coad","_diff",".txt"),row.names=F,col.names=T,sep="\t")
pdf(paste0("all54-coad", "_diff", ".pdf"),width=32,height=32)
forestPlot(mafCompareRes = os01, pVal=0.01, geneFontSize = 0.7)
dev.off()

#oncodrive
laml.sig = oncodrive(maf = laml.plus.gistic3, AACol = 'aaChange', minMut = 10, pvalMethod = 'zscore')
pdf(paste0("0503all54", "_oncodrive.pdf"),width=8,height=8)
plotOncodrive(res = laml.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 0.5)
write.table(laml.sig, file=paste0("0503all54", "_oncodrive.txt"), sep="\t")
dev.off()
laml.sig2 = oncodrive(maf = coad.maf2, AACol = 'HGVSp_Short', minMut = 50, pvalMethod = 'zscore')
pdf(paste0("Coad", "_oncodrive12.pdf"),width=8,height=8)
plotOncodrive(res = laml.sig2, fdrCutOff = 0.05, useFraction = TRUE, labelSize = 0.5)
write.table(laml.sig2, file=paste0("Coad", "_oncodrive0.01.txt"), sep="\t")
dev.off()
laml.sig3 = oncodrive(maf = laml_os3, AACol = 'aaChange', minMut = 2, pvalMethod = 'zscore')
pdf(paste0("0503pri11", "_oncodrive2.pdf"),width=8,height=8)
plotOncodrive(res = laml.sig3, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 0.5)
write.table(laml.sig3, file=paste0("0503pri11", "_oncodrive0.1.txt"), sep="\t")
dev.off()
laml.sig4 = oncodrive(maf = laml_os4, AACol = 'aaChange', minMut = 2, pvalMethod = 'zscore')
pdf(paste0("0503OM15", "_oncodrive3.pdf"),width=8,height=8)
plotOncodrive(res = laml.sig4, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 0.5)
write.table(laml.sig4, file=paste0("0503OM15", "_oncodrive0.1.txt"), sep="\t")
dev.off()

##pri14-ova18-differ
sam7="0504pri11-ova15"
pdf(paste0(sam7, "_RNF43", ".pdf"),width=8,height=8)
lollipopPlot2(m1 = laml_os3, m2 = laml_os4, gene = "RNF43", AACol1 = "aaChange", AACol2 = "aaChange", m1_name = "Pri11", m2_name = "OvaM15")
dev.off()
pdf(paste0(sam7, "_p53", ".pdf"),width=8,height=8)
lollipopPlot2(m1 = laml_os3, m2 = laml_os4, gene = "TP53", AACol1 = "aaChange", AACol2 = "aaChange", m1_name = "Pri11", m2_name = "OvaM15")
dev.off()
pdf(paste0(sam7, "_APC", ".pdf"),width=8,height=8)
lollipopPlot2(m1 = laml_os3, m2 = laml_os4, gene = "APC", AACol1 = "aaChange", AACol2 = "aaChange", m1_name = "Pri11", m2_name = "OvaM15")
dev.off()
pdf(paste0(sam7, "_KRAS", ".pdf"),width=8,height=8)
lollipopPlot2(m1 = laml_os3, m2 = laml_os4, gene = "KRAS", AACol1 = "aaChange", AACol2 = "aaChange", m1_name = "Pri11", m2_name = "OvaM15")
dev.off()
pdf(paste0(sam7, "_ZNF160", ".pdf"),width=8,height=8)
lollipopPlot2(m1 = laml_os3, m2 = laml_os4, gene = "ZNF160", AACol1 = "aaChange", AACol2 = "aaChange", m1_name = "Pri11", m2_name = "OvaM15")
dev.off()
pdf(paste0(sam7, "_ZNF716", ".pdf"),width=8,height=8)
lollipopPlot2(m1 = laml_os3, m2 = laml_os4, gene = "ZNF716", AACol1 = "aaChange", AACol2 = "aaChange", m1_name = "Pri11", m2_name = "OvaM15")
dev.off()
pdf(paste0(sam7, "_STOML1", ".pdf"),width=8,height=8)
lollipopPlot2(m1 = laml_os3, m2 = laml_os4, gene = "STOML1", AACol1 = "aaChange", AACol2 = "aaChange", m1_name = "Pri11", m2_name = "OvaM15")
dev.off()
pdf(paste0(sam7, "_ACVR2A", ".pdf"),width=8,height=8)
lollipopPlot2(m1 = laml_os3, m2 = laml_os4, gene = "ACVR2A", AACol1 = "aaChange", AACol2 = "aaChange", m1_name = "pri14", m2_name = "ova18")
dev.off()
pdf(paste0(sam7, "_PCDHB12", ".pdf"),width=8,height=8)
lollipopPlot2(m1 = laml_os3, m2 = laml_os4, gene = "PCDHB12", AACol1 = "aaChange", AACol2 = "aaChange", m1_name = "pri14", m2_name = "ova18")
dev.off()
pdf(paste0(sam7, "_NLGN1", ".pdf"),width=8,height=8)
lollipopPlot2(m1 = laml_os3, m2 = laml_os4, gene = "NLGN1", AACol1 = "aaChange", AACol2 = "aaChange", m1_name = "pri14", m2_name = "ova18")
dev.off()
pdf(paste0(sam7, "_SMIM3", ".pdf"),width=8,height=8)
lollipopPlot2(m1 = laml_os3, m2 = laml_os4, gene = "SMIM3", AACol1 = "aaChange", AACol2 = "aaChange", m1_name = "pri14", m2_name = "ova18")
dev.off()
pdf(paste0(sam7, "_BRAF", ".pdf"),width=8,height=8)
lollipopPlot2(m1 = laml_os3, m2 = laml_os4, gene = "BRAF", AACol1 = "aaChange", AACol2 = "aaChange", m1_name = "pri14", m2_name = "ova18")
dev.off()

pdf(paste0(sam, "_muc16", ".pdf"),width=8,height=8)
lollipopPlot2(m1 = laml_os3, m2 = laml_os4, gene = "MUC16", AACol1 = "aaChange", AACol2 = "aaChange", m1_name = m1, m2_name = m2)
dev.off()
pdf(paste0(sam, "_fat4", ".pdf"),width=8,height=8)
lollipopPlot2(m1 = laml_os1, m2 = laml_os0, gene = "FAT4", AACol1 = "aaChange", AACol2 = "aaChange", m1_name = m1, m2_name = m2)
dev.off()
pdf(paste0(sam, "_pppIR3A", ".pdf"),width=8,height=8)
lollipopPlot2(m1 = laml_os1, m2 = laml_os0, gene = "PPP1R3A", AACol1 = "aaChange", AACol2 = "aaChange", m1_name = m1, m2_name = m2)
dev.off()
pdf(paste0(sam, "_DMD", ".pdf"),width=8,height=8)
lollipopPlot2(m1 = laml_os1, m2 = laml_os0, gene = "DMD", AACol1 = "aaChange", AACol2 = "aaChange", m1_name = m1, m2_name = m2)
dev.off()
pdf(paste0(sam, "_MUC12", ".pdf"),width=8,height=8)
lollipopPlot2(m1 = laml_os1, m2 = laml_os0, gene = "MUC12", AACol1 = "aaChange", AACol2 = "aaChange", m1_name = m1, m2_name = m2)
dev.off()
pdf(paste0(sam, "_DYNC1H1", ".pdf"),width=8,height=8)
lollipopPlot2(m1 = laml_os1, m2 = laml_os0, gene = "DYNC1H1", AACol1 = "aaChange", AACol2 = "aaChange", m1_name = m1, m2_name = m2)
dev.off()
pdf(paste0(sam, "_DNAH6", ".pdf"),width=8,height=8)
lollipopPlot2(m1 = laml_os1, m2 = laml_os0, gene = "DNAH6", AACol1 = "aaChange", AACol2 = "aaChange", m1_name = m1, m2_name = m2)
dev.off()
pdf(paste0(sam, "_LRP1B", ".pdf"),width=8,height=8)
lollipopPlot2(m1 = laml_os1, m2 = laml_os0, gene = "LRP1B", AACol1 = "aaChange", AACol2 = "aaChange", m1_name = m1, m2_name = m2)
dev.off()

##pathway
sam5="0504all54"
pdf(paste0(sam5,"_pathway.pdf"),width=8,height=8)
OncogenicPathways(maf = laml.plus.gistic3, fontSize=1.2)
dev.off()
pdf(paste0(sam5,"_RAS.pdf"),width=16,height=8)
PlotOncogenicPathways(maf = laml.plus.gistic3, pathways = "RTK-RAS", showTumorSampleBarcodes=T, removeNonMutated=FALSE)
dev.off()
pdf(paste0(sam5,"_NOTCH.pdf"),width=16,height=8)
PlotOncogenicPathways(maf = laml.plus.gistic3, pathways = "NOTCH", showTumorSampleBarcodes=T, removeNonMutated=FALSE)
dev.off()
pdf(paste0(sam5,"_Hippo.pdf"),width=16,height=10)
PlotOncogenicPathways(maf = laml.plus.gistic3, pathways = "Hippo", showTumorSampleBarcodes=T, removeNonMutated=FALSE)
dev.off()
pdf(paste0(sam5,"_WNT.pdf"),width=16,height=8)
PlotOncogenicPathways(maf = laml.plus.gistic3, pathways = "WNT", showTumorSampleBarcodes=T, removeNonMutated=FALSE)
dev.off()
pdf(paste0(sam5,"_PI3K.pdf"),width=16,height=8)
PlotOncogenicPathways(maf = laml.plus.gistic3, pathways = "PI3K", showTumorSampleBarcodes=T, removeNonMutated=FALSE)
dev.off()

pdf(paste0("0611all54","_NRF2.pdf"),width=16,height=8)
PlotOncogenicPathways(maf = laml.plus.gistic3, pathways = "NRF2", showTumorSampleBarcodes=T, removeNonMutated=F, fontSize=0.8, SampleNamefontSize=0.8)
dev.off()

pdf(paste0("0611all54","_Cell_Cycle_0611.pdf"),width=16,height=8)
PlotOncogenicPathways(maf = laml.plus.gistic3, pathways = "Cell_Cycle", showTumorSampleBarcodes=T, removeNonMutated=F, fontSize=0.8, SampleNamefontSize=0.8)
dev.off()

sam5="0504COAD"
mafpath<-coad.maf2
pdf(paste0(sam5,"_pathway.pdf"),width=8,height=8)
OncogenicPathways(maf = mafpath, fontSize=1.2)
dev.off()
pdf(paste0(sam5,"_RAS.pdf"),width=16,height=8)
PlotOncogenicPathways(maf = mafpath, pathways = "RTK-RAS", showTumorSampleBarcodes=F, removeNonMutated=T)
dev.off()
pdf(paste0(sam5,"_NOTCH.pdf"),width=16,height=8)
PlotOncogenicPathways(maf = mafpath, pathways = "NOTCH", showTumorSampleBarcodes=F, removeNonMutated=T)
dev.off()
pdf(paste0(sam5,"_Hippo.pdf"),width=16,height=10)
PlotOncogenicPathways(maf = mafpath, pathways = "Hippo", showTumorSampleBarcodes=F, removeNonMutated=T)
dev.off()
pdf(paste0(sam5,"_WNT.pdf"),width=16,height=8)
PlotOncogenicPathways(maf = mafpath, pathways = "WNT", showTumorSampleBarcodes=F, removeNonMutated=T)
dev.off()
pdf(paste0(sam5,"_PI3K.pdf"),width=16,height=8)
PlotOncogenicPathways(maf = mafpath, pathways = "PI3K", showTumorSampleBarcodes=F, removeNonMutated=T)
dev.off()

sam5="0504pri11"
mafpath<-laml_os3
pdf(paste0(sam5,"_pathway.pdf"),width=8,height=8)
OncogenicPathways(maf = mafpath, fontSize=1.2)
dev.off()
pdf(paste0(sam5,"_RAS.pdf"),width=16,height=8)
PlotOncogenicPathways(maf = mafpath, pathways = "RTK-RAS", showTumorSampleBarcodes=T, removeNonMutated=F)
dev.off()
pdf(paste0(sam5,"_NOTCH.pdf"),width=16,height=8)
PlotOncogenicPathways(maf = mafpath, pathways = "NOTCH", showTumorSampleBarcodes=T, removeNonMutated=F)
dev.off()
pdf(paste0(sam5,"_Hippo.pdf"),width=16,height=10)
PlotOncogenicPathways(maf = mafpath, pathways = "Hippo", showTumorSampleBarcodes=T, removeNonMutated=F)
dev.off()
pdf(paste0(sam5,"_WNT.pdf"),width=16,height=8)
PlotOncogenicPathways(maf = mafpath, pathways = "WNT", showTumorSampleBarcodes=T, removeNonMutated=F)
dev.off()
pdf(paste0(sam5,"_PI3K.pdf"),width=16,height=8)
PlotOncogenicPathways(maf = mafpath, pathways = "PI3K", showTumorSampleBarcodes=T, removeNonMutated=F)
dev.off()

sam5="0504OvaM15"
mafpath<-laml_os4
pdf(paste0(sam5,"_pathway.pdf"),width=8,height=8)
OncogenicPathways(maf = mafpath, fontSize=1.2)
dev.off()
pdf(paste0(sam5,"_RAS.pdf"),width=16,height=8)
PlotOncogenicPathways(maf = mafpath, pathways = "RTK-RAS", showTumorSampleBarcodes=T, removeNonMutated=F, fontSize=0.8, SampleNamefontSize=0.8)
dev.off()
pdf(paste0(sam5,"_NOTCH.pdf"),width=16,height=8)
PlotOncogenicPathways(maf = mafpath, pathways = "NOTCH", showTumorSampleBarcodes=T, removeNonMutated=F, fontSize=0.8, SampleNamefontSize=0.8)
dev.off()
pdf(paste0(sam5,"_Hippo.pdf"),width=16,height=10)
PlotOncogenicPathways(maf = mafpath, pathways = "Hippo", showTumorSampleBarcodes=T, removeNonMutated=F, fontSize=0.8, SampleNamefontSize=0.8)
dev.off()
pdf(paste0(sam5,"_WNT.pdf"),width=16,height=8)
PlotOncogenicPathways(maf = mafpath, pathways = "WNT", showTumorSampleBarcodes=T, removeNonMutated=F, fontSize=0.8, SampleNamefontSize=0.8)
dev.off()
pdf(paste0(sam5,"_PI3K.pdf"),width=16,height=8)
PlotOncogenicPathways(maf = mafpath, pathways = "PI3K", showTumorSampleBarcodes=T, removeNonMutated=F, fontSize=0.8, SampleNamefontSize=0.8)
dev.off()
pdf(paste0(sam5,"_tgf.pdf"),width=16,height=8)
PlotOncogenicPathways(maf = mafpath, pathways = "TGF-Beta", showTumorSampleBarcodes=T, removeNonMutated=F, fontSize=0.8, SampleNamefontSize=0.8)
dev.off()

pdf(paste0(sam5,"_MYC.pdf"),width=16,height=8)
PlotOncogenicPathways(maf = mafpath, pathways = "MYC", showTumorSampleBarcodes=T, removeNonMutated=F, fontSize=0.8, SampleNamefontSize=0.8)
dev.off()
pdf(paste0(sam5,"_Tp53.pdf"),width=16,height=8)
PlotOncogenicPathways(maf = mafpath, pathways = "TP53", showTumorSampleBarcodes=T, removeNonMutated=F, fontSize=0.8, SampleNamefontSize=0.8)
dev.off()

pdf(paste0(sam5,"_TGF-Beta_0611.pdf"),width=16,height=8)
PlotOncogenicPathways(maf = mafpath, pathways = "TGF-Beta", showTumorSampleBarcodes=T, removeNonMutated=F, fontSize=0.8, SampleNamefontSize=0.8)
dev.off()

pdf(paste0(sam5,"_PI3K_0611.pdf"),width=16,height=8)
PlotOncogenicPathways(maf = mafpath, pathways = "PI3K", showTumorSampleBarcodes=T, removeNonMutated=F, fontSize=0.8, SampleNamefontSize=0.8)
dev.off()



##all54-pri-ova-differ
sam6="all54-pri-ova"
laml.plus.gistic3<-get(load('0502all54-cnv-snv-clinical.Rdata'))
laml_os3 = subsetMaf(maf = laml.plus.gistic3, clinQuery = "metaSite %in% 'Pri'")
laml_os4 = subsetMaf(maf = laml.plus.gistic3, clinQuery = "metaSite %in% 'OM'")
save(laml_os3,file='all54-Primary-cnv-snv-clinical.Rdata')
save(laml_os4,file='all54-OVA-cnv-snv-clinical.Rdata')
os02 <- mafCompare(m1 = laml_os3, m2 = laml_os4, m1Name = 'Primary11', m2Name = 'OvaryM15', minMut = 5, useCNV=TRUE)
print(os02)
write.table(os02$result,file=paste0(sam6,"_diff",".txt"),row.names=F,col.names=T,sep="\t")
pdf(paste0(sam6, "_diff", ".pdf"),width=32,height=32)
forestPlot(mafCompareRes = os02, pVal=0.15, geneFontSize = 0.7)
dev.off()
pdf(paste0("0611Pri_ovary_differ_driver24", ".pdf"),width=8,height=8)
coOncoplot(m1 = laml_os3, m2 = laml_os4, m1Name = 'Primary', m2Name = 'OvaryM', removeNonMutated = F, geneNamefont=0.8, genes=driver_genes24)
dev.off()
pdf(paste0("0611pri_Ova_cobar_driver24", ".pdf"),width=8,height=8)
coBarplot(m1 = laml_os3, m2 = laml_os4, m1Name = 'Primary', m2Name = 'OvaryM', geneSize=0.8, genes=driver_genes24)
dev.off()

##all54-pri-PM-differ
sam7="all54-pri-pm"
laml.plus.gistic3<-get(load('0502all54-cnv-snv-clinical.Rdata'))
laml_os3 = subsetMaf(maf = laml.plus.gistic3, clinQuery = "metaSite %in% 'Pri'")
laml_os5 = subsetMaf(maf = laml.plus.gistic3, clinQuery = "metaSite %in% 'PM'")
#save(laml_os3,file='all54-Primary-cnv-snv-clinical.Rdata')
save(laml_os5,file='all54-PeriM-cnv-snv-clinical.Rdata')
os02 <- mafCompare(m1 = laml_os3, m2 = laml_os5, m1Name = 'Primary11', m2Name = 'PeriM', minMut = 5, useCNV=TRUE)
print(os02)
write.table(os02$result,file=paste0(sam7,"_diff",".txt"),row.names=F,col.names=T,sep="\t")
pdf(paste0(sam7, "_diff", ".pdf"),width=32,height=32)
forestPlot(mafCompareRes = os02, pVal=0.15, geneFontSize = 0.7)
dev.off()
pdf(paste0("0611Pri_PM_differ_driver21", ".pdf"),width=8,height=8)
coOncoplot(m1 = laml_os3, m2 = laml_os5, m1Name = 'Primary', m2Name = 'PeriM', removeNonMutated = F, geneNamefont=0.8, genes=driver_genes21)
dev.off()
pdf(paste0("0611pri_PM_differ_driver21_coBAR2", ".pdf"),width=8,height=8)
coBarplot(m1 = laml_os3, m2 = laml_os5, m1Name = 'Primary', m2Name = 'PeriM', geneSize=0.8, genes=driver_genes21)
dev.off()

pdf(paste0("0611ova_PM_differ_driver22", ".pdf"),width=8,height=8)
coOncoplot(m1 = laml_os4, m2 = laml_os5, m1Name = 'OvaM15', m2Name = 'PeriM8', removeNonMutated = F, geneNamefont=0.8, genes=driver_genes22)
dev.off()
pdf(paste0("0611ova_PM_cobar_driver22", ".pdf"),width=8,height=8)
coBarplot(m1 = laml_os4, m2 = laml_os5, m1Name = 'OvaM15', m2Name = 'PeriM', geneSize=0.7, genes=driver_genes22)
dev.off()

laml_os6 = subsetMaf(maf = laml.plus.gistic3, clinQuery = "metaSite %in% 'LM'")
#save(laml_os6,file='all54-lymphM-cnv-snv-clinical.Rdata')
pdf(paste0("0611ova_lymphM_differ_driver22", ".pdf"),width=8,height=8)
coOncoplot(m1 = laml_os4, m2 = laml_os6, m1Name = 'OvaM15', m2Name = 'LymphM', removeNonMutated = F, geneNamefont=0.8, genes=driver_genes22)
dev.off()
pdf(paste0("0611ova_lymphM_cobar_driver22", ".pdf"),width=8,height=8)
coBarplot(m1 = laml_os4, m2 = laml_os6, m1Name = 'OvaM15', m2Name = 'LymphM', geneSize=0.7, genes=driver_genes22)
dev.off()

laml_os7 = subsetMaf(maf = laml.plus.gistic3, clinQuery = "metaSite %in% 'OME'")
#save(laml_os6,file='all54-Omentum-cnv-snv-clinical.Rdata')
pdf(paste0("0611ova_omentumM_differ_driver24", ".pdf"),width=8,height=8)
coOncoplot(m1 = laml_os4, m2 = laml_os7, m1Name = 'OvaM15', m2Name = 'OmenM', removeNonMutated = F, geneNamefont=0.6, genes=driver_genes24)
dev.off()
pdf(paste0("0611ova_omentumM_differ_driver24", ".pdf"),width=8,height=8)
coBarplot(m1 = laml_os4, m2 = laml_os7, m1Name = 'OvaM15', m2Name = 'OmenM', geneSize=0.6, genes=driver_genes24)
dev.off()

pdf(paste0("0611_all54_coad_differ_drivergenes22", ".pdf"),width=16,height=8)
coOncoplot(m1 = laml.plus.gistic2, m2 = coad.maf2, m1Name = 'All54', m2Name = 'COAD', removeNonMutated = T, geneNamefont=0.6, genes=driver_genes22)
dev.off()
pdf(paste0("0611_all54_coad_differ_drivergenes22_coBar", ".pdf"),width=8,height=8)
coBarplot(m1 = laml.plus.gistic2, m2 = coad.maf2, m1Name = 'All54', m2Name = 'COAD', geneSize=0.6, genes=driver_genes22)
dev.off()

##survival
setwd('0330all54-maf')
coad.maf2<-get(load('TCGA-coad-cnv-snv-clinical.Rdata'))
coadall.lesions="202203TCGA-coad-gistic/all_lesions.conf_99.txt"
coadamp.genes = "202203TCGA-coad-gistic/amp_genes.conf_99.txt"
coaddel.genes = "202203TCGA-coad-gistic/del_genes.conf_99.txt"
coadscores.gis = "202203TCGA-coad-gistic/scores.gistic"
clin <- read.csv("./202203TCGA-coad-gistic/TCGA-COAD_Clinical.csv", header=T)

#TCGA-snp-clinical
coad<-TCGAmutations::tcga_load(study = "COAD")
save(coad,file='TCGA-coad-snv.Rdata')
coad2<-get(load('TCGA-coad-snv.Rdata'))
head(coad2@data)

library(TCGAbiolinks)
query1 <- GDCquery( project = "TCGA-COAD", data.category = "Simple Nucleotide Variation", data.type = "Masked Somatic Mutation", legacy=F)
#GDCdownload(query1, directory = "GDCdata/")
coadmuts <- GDCprepare(query1, directory = "GDCdata/")
clin<-GDCquery_clinic(project='TCGA-COAD', type='clinical')
head(clin)
write.csv(clin,file = "TCGA-COAD-clin.csv")
clin<-read.csv(file = "TCGA-COAD-clin.csv")

pdf(paste0("0912TCGA-COAD", "_snv_summary_top50.pdf"),width=8,height=8)
plotmafSummary(maf = coad.maf2, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE, showBarcodes=F, top=50)
dev.off()

clin2 <- clin %>%
  dplyr::select(c(
    "submitter_id", "days_to_last_follow_up", "vital_status",
    "ajcc_pathologic_t", "ajcc_pathologic_n", "ajcc_pathologic_m",
    "ajcc_pathologic_stage", "gender"
    )) %>%
  `names<-`(c("Tumor_Sample_Barcode", "time", "status", "T", "N", "M", "stage", "gender"))
#coad.maf <- read.maf(maf=coadmuts,
  #clinicalData = clin2,
  #isTCGA = TRUE)

clin2[is.na(clin2)]=0
clin2$status<-ifelse(clin2$status == "Alive",1,0)
#colnames(clin)[1]  = 'Tumor_Sample_Barcode'
head(clin)

coad.maf3 = read.maf(
maf = coadmuts,
  gisticAllLesionsFile = coadall.lesions,
  gisticAmpGenesFile = coadamp.genes,
  gisticDelGenesFile = coaddel.genes,
  gisticScoresFile = coadscores.gis,
  isTCGA = T,
  verbose = FALSE,
  clinicalData = clin2
)
getClinicalData(coad.maf3)

coad.maf3@data$Tumor_Sample_Barcode<-substr(coad.maf2@data$Tumor_Sample_Barcode,1,12)
head(coad.maf3@data$Tumor_Sample_Barcode)

laml.plus.gistic3<-get(load("0502all54-cnv-snv-clinical.Rdata"))
getClinicalData(laml.plus.gistic3)
coad.maf3<-get(load("TCGA-coad-cnv-snv-clin-0618.Rdata"))
getClinicalData(coad.maf3)

coad.maf4 = read.maf(
+ maf = coad.maf3@data,
+     isTCGA = T,
+   verbose = FALSE,
+   clinicalData = clin2
+ )
save(coad.maf4,file="0620coad_snv_CNV_clin.Rdata")


laml_os3 = subsetMaf(maf = laml.plus.gistic3, clinQuery = "metaSite %in% 'Pri'")
laml_os4 = subsetMaf(maf = laml.plus.gistic3, clinQuery = "metaSite %in% 'OM'")
laml_os1 = subsetMaf(maf = coad.maf4, clinQuery = "M %in% 'M0'")
#laml_os0 = subsetMaf(maf = coad.maf2, clinQuery = "!(M %in% 'M0')")
laml_os0 = subsetMaf(maf = coad.maf4, clinQuery = "M %in% c('M1','M1a')")

###Predict genesets associated with survival-COAD_M1a
coad_m1<-survGroup(
maf = coad.maf4,
  top = 30,
  genes = NULL,
  geneSetSize = 3,
  minSamples = 3,
  clinicalData = NULL,
  time = "time",
  Status = "status",
  verbose = TRUE,
  plot = F
)
write.csv(coad_m1,file = "coad-survival-genes2.csv")

Geneset:  APC,TP53,KRAS,TTN [N= 8 ]
Geneset:  APC,TP53,KRAS,PIK3CA [N= 5 ]
       Gene_combination P_value       hr WT Mutant
1:    APC_TP53_KRAS_TTN  0.0455 1.04e-08 54      8
2: APC_TP53_KRAS_PIK3CA  0.0899 2.04e-01 57      5
 Gene_combination P_value       hr WT Mutant
1:  APC_TP53_ADGRB3  0.0240 5.02e+00 59      3
2:   APC_TP53_CSMD2  0.0414 1.51e-01 58      4
3:    TP53_KRAS_TTN  0.0455 1.04e-08 54      8
4:  APC_SYNE1_CSMD3  0.0463 3.38e+00 58      4
5:     APC_TP53_CFH  0.0540 7.37e+00 59      3

1:        APC_PREX2 0.000277 17.900 58      4
2:        APC_CHRM2 0.006310  4.440 58      4
3:        APC_CSMD2 0.041400  0.151 58      4
4:        APC_CSMD3 0.066200  2.340 54      8
5:       TP53_CSMD2 0.076400  0.282 57      5

pdf(paste0("0620-coad-M0-survival-ADGRB3", ".pdf"),width=8,height=8)
mafSurvGroup(maf = laml_os1, geneSet = c("ADGRB3", "TP53", "APC"), time = 'time', Status = 'status')
dev.off()
 
#DST_BRAF_RNF43
pdf(paste0("0620-coad-survival-DST_BRAF_RNF43", ".pdf"),width=8,height=8)
mafSurvGroup(maf = laml_os1, geneSet = c("DST", "BRAF", "RNF43"), time = 'time', Status = 'status')
dev.off()

pdf(paste0("0620-coad-M1a-survival-TP53-RNF43", ".pdf"),width=8,height=8)
mafSurvGroup(maf = laml_os0, geneSet = c("RNF43", "TP53"), time = 'time', Status = 'status')
dev.off()

driver_genes3= c("TP53", "BSN", "RYR3")
mafSurvGroup(maf = laml.plus.gistic3, geneSet = driver_genes3, time = 'days_to_last_followup', Status = 'OS12')

driver_genes3= c("TP53", "NRAS")
mafSurvGroup(maf = laml_os3, geneSet = driver_genes3, time = 'days_to_last_followup', Status = 'OS24')
Looking for clinical data in annoatation slot of MAF..
Median survival..
    Group medianTime N
1: Mutant        270 5
2:     WT        540 6
 
driver_genes3= c("TP53", "NRAS", "DVL2")
 mafSurvGroup(maf = laml_os3, geneSet = driver_genes3, time = 'days_to_last_followup', Status = 'OS24')
Looking for clinical data in annoatation slot of MAF..
Median survival..
    Group medianTime N
1: Mutant        300 4
2:     WT        360 7

mafSurvival(maf = coad.maf3, genes = "NRAS", time = 'time', Status = 'status', isTCGA = TRUE)

###Predict genesets associated with survival
sur_54<-survGroup(
  maf = laml.plus.gistic3,
  top = 20,
  genes = NULL,
  geneSetSize = 2,
  minSamples = 5,
  clinicalData = NULL,
  time = "days_to_last_followup",
  Status = "OS24",
  verbose = TRUE,
  plot = T
)
write.csv(sur_54,file = "all54-surv5-genes2.csv")
##ova
sur_ova<-survGroup(
   maf = laml_os4,
   top = 100,
   genes = NULL,
   geneSetSize = 2,
  minSamples = 2,
   clinicalData = NULL,
   time = "days_to_last_followup",
   Status = "OS12",
   verbose = TRUE,
   plot = F
 )
write.csv(sur_ova,file = "ova15-surv-genes2.csv")

mafSurvGroup(maf = laml_os4, geneSet = "RNF43", time = 'days_to_last_followup', Status = 'OS24')
mafSurvival(maf = laml_os4, genes = "RNF43", time = 'days_to_last_followup', Status = 'OS24')
dev.off()

###OVA-Metastasis-gene2-survival
#genes=c("MUC16","ARID1A","PRKDC")
pdf(paste0("coad-460-survival0704-", "MUC16_ARID1A_PRKDC", ".pdf"),width=8,height=8)
mafSurvGroup(maf = coad.maf4, geneSet = genes, time = 'time', Status = 'status')
genes=c("TP53","ARID1A","PRKDC")
pdf(paste0("0701-OVA-survival2-", genes, ".pdf"),width=8,height=8)
mafSurvGroup(maf = laml_os4, geneSet = genes, time = 'days_to_last_followup', Status = 'OS12')
dev.off()
pdf(paste0("0701-all54-pri11-survival3-", genes, ".pdf"),width=8,height=8)
mafSurvGroup(maf = laml_os3, geneSet = genes, time = 'days_to_last_followup', Status = 'OS12')
dev.off()
pdf(paste0("0701-coad-460-survival3-", genes, ".pdf"),width=8,height=8)
mafSurvGroup(maf = coad.maf4, geneSet = genes, time = 'time', Status = 'status')
dev.off()
pdf(paste0("0701-all54-survival3-", genes, ".pdf"),width=8,height=8)
mafSurvGroup(maf = laml.plus.gistic3, geneSet = genes, time = 'days_to_last_followup', Status = 'OS12')
dev.off()

