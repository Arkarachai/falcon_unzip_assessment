library(ggplot2)
library(stats)

args = commandArgs(trailingOnly=TRUE)
dataset=read.table(args[1])
high=dataset[dataset$V1=='mis_join',]
low=dataset[dataset$V1=='random_error',]
colnames(dataset)=c('label','value')
boxplot(high$V2,low$V2)
boxplot(high$V2,low$V2,log='y')
dataset$label=as.factor(dataset$label)

pdf(paste(args[1],'.pdf'))
ggplot(dataset,aes(x=label,y=value),log='y')+geom_violin()+geom_boxplot(width=.12, fill=I('black'), notch=T, col='grey40')+ylab('SNPs distance at boundary region')+scale_x_discrete(labels=c("Misjoined haplotig","Random error haplotig"))+xlab(NULL)+coord_cartesian(ylim = c(0, 60000)) 
dev.off()

#t.test(dataset$value~dataset$label) # normal distribute only; maybe I should try permutation test or bootstrap test
wilcox.test(value~label,data=dataset)
summary(high$V2)
summary(low$V2)
IQR(high$V2)
IQR(low$V2)
ansari.test(high$V2, low$V2)
var(high$V2)
var(low$V2)
var.test(high$V2,low$V2) 
