##gisitic2maftools
R
setwd('/mydata/202202_all54_allseg2')
rm(list = ls())
options(stringsAsFactors = F)
library(maftools)
laml.gistic = readGistic(gisticAllLesionsFile = 'all_lesions.conf_90.txt',
                         gisticAmpGenesFile = 'amp_genes.conf_90.txt', 
                         gisticDelGenesFile = 'del_genes.conf_90.txt', 
                         gisticScoresFile = 'scores.gistic', 
                         isTCGA = T)
jpeg(file = "female54_maftools.jpeg")
pdf(file = "female54_maftools.pdf")
gisticChromPlot(gistic = laml.gistic, ref.build = 'hg38')
dev.off()
 
jpeg(file = "feamle54_maftools_bubble.jpeg")
gisticBubblePlot(gistic = laml.gistic)
dev.off()

jpeg(file = "female54_oncoplot_maftools.jpeg")
gisticOncoPlot(gistic = laml.gistic, sortByAnnotation = TRUE, top = 30)
dev.off()
 
