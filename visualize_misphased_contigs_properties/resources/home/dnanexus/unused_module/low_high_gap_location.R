args = commandArgs(trailingOnly=TRUE)
dataset=read.table(args[1])
high=dataset[dataset$V1=='mis_join',]
low=dataset[dataset$V1=='random_error',]

pdf(paste(args[1],'.pdf'))
par(mfrow=c(2,1))
hist(high$V2,xlab=NULL,freq=F,breaks=20,main="Mis-joined haplotigs")
hist(low$V2,xlab='Normalized coordinates',freq=F,breaks=20,main="Random error haplotigs")
dev.off()