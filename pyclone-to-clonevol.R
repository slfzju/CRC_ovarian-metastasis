
##MACPRO_terminal
#! /path/to/Rscript
##pyclon-tsv-merg-clonevol-tsv
ls *_loci.tsv>config1;
cat config1 | while read id;
do
        i=${id};
         echo $i
cat $i | cut -f 1-6 | sed '1d' | paste - - - - > ${i:0:9}_loci2.tsv;
cat *loci2.tsv > gaozh_all_loci3.tsv;
done

R --max-ppsize 500000
install.packages('data.table', repos = 'https://mirrors.tuna.tsinghua.edu.cn/CRAN')
setwd("20221012lal_tree")

options("expressions"=500000)
options("expression"=500000)
rm(list=ls())
gc()
gc(verbose=TRUE, full=TRUE)
library(clonevol)
library(data.table)
library(clonevol)
library(tidyverse)
library(ggplot2)
getwd()

yuli<-as.data.frame(fread("cjl_all1220_clu16.tsv", sep="\t", header=F, dec=".", fill=TRUE))
names(yuli)<-c('mutation', 'sample', 'cluster', 'cellular_prevalence', 'cellular_prevalence_std', 'variant_allele_frequency')
yuli$vaf <-yuli$cellular_prevalence*100/2
yuli<-as.data.table(yuli)
yuli3 <- dcast(yuli, formula = mutation + cluster ~ sample, fun.aggregate = mean, value.var = "vaf")
yuli3[is.na(yuli3)] <- 0
yuli3 <- yuli3[order(yuli3$cluster),]
yuli3$cluster
data2<-data.frame(yuli3)

yuli3$genes<-substring(yuli3$mutation, 15)
data<-data.frame(yuli3)
data$genes<-gsub('^[1-9]','', data$gene)
data$genes<- sub(':', '', data$gene)

##driver_genes
mutsig<-read.table("results.sig_genes.txt", sep = "\t", header = T,nrow=0)
mutsig5<-mutsig[mutsig$p < 0.1,]
mutsig6<-mutsig5[,1]
mutsig6=unique(mutsig6)
mutsig6=unique(mutsig6)
i<-mutsig6
data2 = data %>% filter(data$genes %in% i)
head(data2)
data2$driver <- as.logical(FALSE)
j<-mutsig6
for (j in mutsig6[1:4]){
data2$driver[data2$gene==j] <-as.logical(TRUE)	
return
}
head(data2)
data2$cluster

data2$cluster <-data2$cluster+1
data2 <- data2[order(data2$cluster),]
data2$cluster
data2[data2$cluster == 5,2] <-2
data2[data2$cluster == 7,2] <-6
data2[data2$cluster == 8,2] <-7
write.table(data2,file = "cjl4_all1222_clu.txt", quote= F, sep = "\t", row.names = F)
data4<-read.table("cjl5_vaf1222_clu6_10.txt", sep = "\t", header = T)

colnames(data4)<-c("mutation","cluster","Pri.vaf", "Lym.vaf", "Ova.vaf","Per.vaf")
vaf.col.names <- grep('.vaf', colnames(data4), value=T)
sample.names <- gsub('.vaf', '', vaf.col.names)
data4[, sample.names] <- data4[, vaf.col.names]
vaf.col.names <- sample.names
sample.groups <- c("Pri", "Lym", "Ova", "Per")
names(sample.groups) <- vaf.col.names
sample.groups
data4 <- data4[order(data4$cluster),]
data4<-data.frame(data4)
clone.colors <- NULL

y <- infer.clonal.models(variants = data4, 
cluster.col.name = 'cluster', 
vaf.col.names = vaf.col.names, 
sample.groups = sample.groups, 
cancer.initiation.model = 'polyclonal', 
subclonal.test = 'bootstrap', 
subclonal.test.model = 'non-parametric', 
num.boots = 1000, 
founding.cluster = 1, 
clone.colors = clone.colors,
ignore.clusters = NULL, 
sum.p = 0.05,
alpha = 0.05,
min.cluster.vaf = 0.01)

data4[data4$driver,],
cluster.col.name = data4$cluster,
event.col.name = data4$gene)


y <- convert.consensus.tree.clone.to.branch(y, branch.scale = 'sqrt')
sam='yuli5_1223_1'
outdir= paste0(sam, '_tree')
pdf(file=paste0(sam, '_trees.pdf'), width = 9, height = 15, useDingbats = FALSE)
plot.all.trees.clone.as.branch(y_1, max.num.models.to.plot = 500,
#text.size = 0.8,
branch.width = 0.5,
#branch.angle = 80, 
node.size = 2,
node.label.size = 1,
node.text.size = 1) 
plot.clonal.models(y,  
max.num.models.to.plot=500,
  # box plot parameters
box.plot = TRUE,
fancy.boxplot = TRUE,
fancy.variant.boxplot.highlight.shape = 21,
fancy.variant.boxplot.highlight.fill.color = 'red',
fancy.variant.boxplot.highlight.color = 'black',
fancy.variant.boxplot.highlight.note.col.name = 'gene',
fancy.variant.boxplot.highlight.note.color = 'blue',
fancy.variant.boxplot.highlight.note.size = 2,
fancy.variant.boxplot.jitter.alpha = 1,
fancy.variant.boxplot.jitter.center.color = 'grey50',
fancy.variant.boxplot.base_size = 12,
fancy.variant.boxplot.plot.margin = 1,
fancy.variant.boxplot.vaf.suffix = '.VAF',
  # bell plot parameters
clone.shape = 'bell',
bell.event = TRUE,
bell.event.label.color = 'blue',
bell.event.label.angle = 60,
clone.time.step.scale = 0.5,
bell.curve.step = 1,
  # node-based consensus tree parameters
merged.tree.plot = TRUE,
tree.node.label.split.character = NULL,
tree.node.shape = 'circle',
tree.node.size = 20,
tree.node.text.size = 1,
merged.tree.node.size.scale = 1,
merged.tree.node.text.size.scale = 1,
merged.tree.cell.frac.ci = FALSE,
  # branch-based consensus tree parameters
merged.tree.clone.as.branch = TRUE,
mtcab.event.sep.char = ',',
mtcab.branch.text.size = 0.8,
mtcab.branch.width = 0.4,
mtcab.branch.angle = 45,
mtcab.node.size = 2,
mtcab.node.label.size = 1,
mtcab.node.text.size = 1,
  # cellular population parameters
cell.plot = TRUE,
num.cells = 100,
cell.border.size = 0.25,
cell.border.color = 'black',
clone.grouping = 'horizontal',
  #meta-parameters
scale.monoclonal.cell.frac = TRUE,
show.score = FALSE,
cell.frac.ci = TRUE,
disable.cell.frac = FALSE,
  # output figure parameters
out.dir = outdir,
out.format = 'pdf',
overwrite.output = TRUE,
width = 18,
height = 8,     
 # vector of width scales for each panel from left to right
panel.widths = c(6,8,6,8,8))
dev.off()

