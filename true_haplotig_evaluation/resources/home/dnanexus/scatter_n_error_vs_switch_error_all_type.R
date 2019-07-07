dataset=read.table('haplotig_evaluation_table',sep='\t',header=T)
library(scales)
pchvalue=16 #16
alphavalue=0.1

lowset=dataset[dataset$phasing_type=='random_error',]
highset=dataset[dataset$phasing_type=='mis_join',]
mediumset=data.frame(rbind(dataset[dataset$phasing_type=='inbetween',],dataset[dataset$phasing_type=='mix',]))
mixset=dataset[dataset$phasing_type=='mix',]
inbetweenset=dataset[dataset$phasing_type=='inbetween',]
scaling_factor=3000 #min(dataset$snp_span_length)

# pdf(paste('scatter_n_error_vs_switch_error_all_type','.pdf',sep=""))
# 
# highset_premelt=data.frame(cbind(highset$error_count,highset$parent_switch_n))
# colnames(highset_premelt)=c("error_count","parent_switch_n")
# highset_aggregate=aggregate(list(numdup=rep(1,nrow(highset_premelt))), highset_premelt, length)
# 
# lowset_premelt=data.frame(cbind(lowset$error_count,lowset$parent_switch_n))
# colnames(lowset_premelt)=c("error_count","parent_switch_n")
# lowset_aggregate=aggregate(list(numdup=rep(1,nrow(lowset_premelt))), lowset_premelt, length)
# 
# mediumset_premelt=data.frame(cbind(mediumset$error_count,mediumset$parent_switch_n))
# colnames(mediumset_premelt)=c("error_count","parent_switch_n")
# mediumset_aggregate=aggregate(list(numdup=rep(1,nrow(mediumset_premelt))), mediumset_premelt, length)
# 
# #plot(highset_aggregate$error_count,highset_aggregate$parent_switch_n,xlab='Number of mis-phased SNPs',ylab='Number of switch error',pch=15,col=alpha('blue',alphavalue),type='p',xlim=c(0,100),ylim=c(0,50),cex=log(highset_aggregate$numdup+1))
# plot(highset_aggregate$error_count,highset_aggregate$parent_switch_n,xlab='Number of mis-phased SNPs',ylab='Number of switch error',pch=15,col=alpha('blue',alphavalue),type='p',xlim=c(0,max(dataset$error_count)),ylim=c(0,max(dataset$parent_switch_n)),cex=log10(highset_aggregate$numdup+1))
# lines(lowset_aggregate$error_count,lowset_aggregate$parent_switch_n,pch=16,col=alpha('red',alphavalue),type='p',cex=log(lowset_aggregate$numdup+1))
# lines(mediumset_aggregate$error_count,mediumset_aggregate$parent_switch_n,pch=15,col=alpha('blue',alphavalue),type='p',cex=log(mediumset_aggregate$numdup+1))
# abline(0,2)
# dev.off()
# 
# pdf(paste('scatter_error_density_vs_switch_error_all_type','.pdf',sep=""))
# 
# highset_premelt=data.frame(cbind(highset$error_count/highset$snp_span_length,highset$parent_switch_n))
# colnames(highset_premelt)=c("error_count","parent_switch_n")
# highset_aggregate=aggregate(list(numdup=rep(1,nrow(highset_premelt))), highset_premelt, length)
# 
# lowset_premelt=data.frame(cbind(lowset$error_count/lowset$snp_span_length,lowset$parent_switch_n))
# colnames(lowset_premelt)=c("error_count","parent_switch_n")
# lowset_aggregate=aggregate(list(numdup=rep(1,nrow(lowset_premelt))), lowset_premelt, length)
# 
# mediumset_premelt=data.frame(cbind(mediumset$error_count/mediumset$snp_span_length,mediumset$parent_switch_n))
# colnames(mediumset_premelt)=c("error_count","parent_switch_n")
# mediumset_aggregate=aggregate(list(numdup=rep(1,nrow(mediumset_premelt))), mediumset_premelt, length)
# 
# #plot(highset_aggregate$error_count,highset_aggregate$parent_switch_n,xlab='mis-phased SNPs density',ylab='Number of switch error',pch=15,col=alpha('blue',alphavalue),type='p',xlim=c(0,max(dataset$error_count/dataset$snp_span_length)),ylim=c(0,max(dataset$parent_switch_n)),cex=log(highset_aggregate$numdup+1))
# #lines(lowset_aggregate$error_count,lowset_aggregate$parent_switch_n,pch=16,col=alpha('red',alphavalue),type='p',cex=log(lowset_aggregate$numdup+1))
# #lines(mediumset_aggregate$error_count,mediumset_aggregate$parent_switch_n,pch=15,col=alpha('blue',alphavalue),type='p',cex=log(mediumset_aggregate$numdup+1))
# plot(dataset$error_count/dataset$snp_span_length,dataset$parent_switch_n,xlab='mis-phased SNPs density',ylab='Number of switch error',pch=16,col=alpha('blue',alphavalue),type='p',xlim=c(0,max(dataset$error_count/dataset$snp_span_length)),ylim=c(0,max(dataset$parent_switch_n)),cex=log10(dataset$snp_span_length/scaling_factor)+0.1)
# abline(h=0)
# abline(v=0)
# 
# dev.off()


pdf(paste('scatter_error_ratio_vs_switch_error_all_type','.pdf',sep=""))

# highset_premelt=data.frame(cbind(highset$error_count/highset$N,highset$parent_switch_n))
# colnames(highset_premelt)=c("error_count","parent_switch_n")
# highset_aggregate=aggregate(list(numdup=rep(1,nrow(highset_premelt))), highset_premelt, length)
# 
# lowset_premelt=data.frame(cbind(lowset$error_count/lowset$N,lowset$parent_switch_n))
# colnames(lowset_premelt)=c("error_count","parent_switch_n")
# lowset_aggregate=aggregate(list(numdup=rep(1,nrow(lowset_premelt))), lowset_premelt, length)
# 
# mediumset_premelt=data.frame(cbind(mediumset$error_count/mediumset$N,mediumset$parent_switch_n))
# colnames(mediumset_premelt)=c("error_count","parent_switch_n")
# mediumset_aggregate=aggregate(list(numdup=rep(1,nrow(mediumset_premelt))), mediumset_premelt, length)

#plot(highset_aggregate$error_count,highset_aggregate$parent_switch_n,xlab='mis-phased SNPs ratio',ylab='Number of switch error',pch=15,col=alpha('blue',alphavalue),type='p',xlim=c(0,max(dataset$error_count/dataset$N)),ylim=c(0,max(dataset$parent_switch_n)),cex=log(highset_aggregate$numdup+1))
#lines(lowset_aggregate$error_count,lowset_aggregate$parent_switch_n,pch=16,col=alpha('red',alphavalue),type='p',cex=log(lowset_aggregate$numdup+1))
#lines(mediumset_aggregate$error_count,mediumset_aggregate$parent_switch_n,pch=15,col=alpha('blue',alphavalue),type='p',cex=log(mediumset_aggregate$numdup+1))
#abline(0,2)
#max(dataset$error_count/dataset$N)
plot(dataset$error_count/dataset$N,dataset$parent_switch_n,xlab='mis-phased SNVs ratio',ylab='Number of switch error',cex.lab=1.25,pch=16,col=alpha('blue',alphavalue),type='p',xlim=c(0,1),ylim=c(0,max(dataset$parent_switch_n)),cex=log10(dataset$snp_span_length/scaling_factor)+0.1)
abline(h=0)
abline(v=0)

dev.off()
