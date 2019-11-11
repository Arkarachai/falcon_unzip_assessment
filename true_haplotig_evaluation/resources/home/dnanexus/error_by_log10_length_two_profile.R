args = commandArgs(trailingOnly=TRUE)
dataset1=read.table('haplotig_evaluation_table1',sep='\t',header=T)
dataset2=read.table('haplotig_evaluation_table2',sep='\t',header=T)

library(scales)

pchvalue=3 #16 
plot_color='blue'
alphavalue=0.3
contig_size_cutoff=1000000

#pdf(paste(args,'_scatter_error_by_log10_length','.pdf',sep=""))
plot(dataset1$error_count/dataset1$N,log10(dataset1$snp_span_length),
     xlab='Mis-phased SNVs Fraction',ylab='Haplotig Size (log10 bp)',
     cex.lab=1.5,cex.axis=1.5,cex=0.75,pch=pchvalue,col=alpha(plot_color,alphavalue),
     type='p',xlim=c(0,1),ylim=c(3,8))
abline(h=log10(contig_size_cutoff),lty=2)

pchvalue=4 #17
plot_color='red'
alphavalue=0.3
contig_size_cutoff=1000000

lines(dataset2$error_count/dataset2$N,log10(dataset2$snp_span_length),
     xlab='Mis-phased SNVs Fraction',ylab='Haplotig Size (log10 bp)',
     cex.lab=1.5,cex.axis=1.5,cex=0.75,pch=pchvalue,col=alpha(plot_color,alphavalue),
     type='p',xlim=c(0,1),ylim=c(3,8))
abline(h=log10(contig_size_cutoff),lty=2)
dev.off()
