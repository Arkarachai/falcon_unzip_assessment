args = commandArgs(trailingOnly=TRUE)
dataset=read.table(args,sep='\t',header=T)
library(scales)

pchvalue=16 
plot_color='blue'
alphavalue=0.3
contig_size_cutoff=1000000

pdf(paste(args,'_scatter_error_by_log10_length','.pdf',sep=""))
plot(dataset$error_count/dataset$N,log10(dataset$snp_span_length),
     xlab='Mis-phased SNVs Fraction',ylab='Haplotig Size (log10 bp)',
     cex.lab=1.5,cex.axis=1.5,cex=0.75,pch=16,col=alpha(plot_color,alphavalue),
     type='p',xlim=c(0,1),ylim=c(3,8))
abline(h=log10(contig_size_cutoff),lty=2)
dev.off()
