R
library('sequenza')
out_dir="/mydata"
sam="all54"
seq_file <- $sam_filter.vcf
#seqz.data <- VarScan2seqz(varscan.somatic = snp)
test <- sequenza.extract(seq_file, chromosome.list=paste("chr",c((1:22),"X","Y"),sep=""), verbose = FALSE)
CP <- sequenza.fit(test)
sequenza.results(sequenza.extract = test, cp.table = CP, sample.id = sam, out.dir=out_dir)
cint <- get.ci(CP)
cellularity <- cint$max.cellularity
write.table(cellularity, paste(out_dir,paste(sam,"_cellularity.txt",sep=""),sep='/'), col.names = FALSE, row.names = FALSE, sep = "\t", quote=FALSE)
ploidy <- cint$max.ploidy
write.table(ploidy, paste(out_dir,paste(sam,"_ploidy.txt",sep=""),sep='/'), col.names = FALSE, row.names = FALSE, sep = "\t", quote=FALSE)
