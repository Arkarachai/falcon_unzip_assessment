dataset=read.table('haplotig_evaluation_table',sep='\t',header=T)
library(scales)
lowset=dataset[dataset$phasing_type=='random_error',]
highset=dataset[dataset$phasing_type=='mis_join',]
mediumset=data.frame(rbind(dataset[dataset$phasing_type=='inbetween',],dataset[dataset$phasing_type=='mix',]))
mixset=dataset[dataset$phasing_type=='mix',]
inbetweenset=dataset[dataset$phasing_type=='inbetween',]
pdf(paste('scatter_n_error_vs_switch_error','.pdf',sep=""))
pchvalue=16 #16
alphavalue=0.2
highset_premelt=data.frame(cbind(highset$error_count,highset$parent_switch_n))
colnames(highset_premelt)=c("error_count","parent_switch_n")
highset_aggregate=aggregate(list(numdup=rep(1,nrow(highset_premelt))), highset_premelt, length)

lowset_premelt=data.frame(cbind(lowset$error_count,lowset$parent_switch_n))
colnames(lowset_premelt)=c("error_count","parent_switch_n")
lowset_aggregate=aggregate(list(numdup=rep(1,nrow(lowset_premelt))), lowset_premelt, length)

plot(highset_aggregate$error_count,highset_aggregate$parent_switch_n,xlab='Number of mis-phased SNPs',ylab='Number of switch error',pch=15,col=alpha('blue',alphavalue),type='p',xlim=c(0,528),ylim=c(0,10),cex=log(highset_aggregate$numdup+1))
lines(lowset_aggregate$error_count,lowset_aggregate$parent_switch_n,pch=16,col=alpha('red',alphavalue),type='p',cex=log(lowset_aggregate$numdup+1))
dev.off()
