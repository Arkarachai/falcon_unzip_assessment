library(ggplot2)

dataset=read.table('boundary_size_vs_location',header=T)
head(dataset)
random_error=dataset[dataset$class=='random_error',]
mis_join=dataset[dataset$class=='mis_join',]

# location
pdf(paste('histogram_normalize_boundary_location','.pdf',sep=""))
par(mfrow=c(2,1))
hist(mis_join$n_end_position,xlab=NULL,freq=F,breaks=20,main="Mis-joined haplotigs")
hist(random_error$n_end_position,xlab='Normalized coordinates',freq=F,breaks=20,main="Random error haplotigs")
dev.off()
hist(mis_join$n_start_position,xlab=NULL,freq=F,breaks=20,main="Mis-joined haplotigs")
hist(random_error$n_start_position,xlab='Normalized coordinates',freq=F,breaks=20,main="Random error haplotigs")

# boundary size
paste('mean misjoin boundary distance',median(mis_join$distance))
paste('mean random_error boundary distance',median(random_error$distance))
#summary(mis_join)
#summary(random_error)
print('difference of median test for misjoin vs random_error')
wilcox.test(distance~class,data=dataset)
paste('IQR misjoin boundary distance',IQR(mis_join$distance))
paste('IQR random error boundary distance',IQR(random_error$distance))
print('equality of variance test for misjoin vs random_error')
ansari.test(mis_join$distance, random_error$distance)
#var.test(mis_join$distance, random_error$distance)
pdf(paste('box_plot_boundary_size','.pdf',sep=""))
ggplot(dataset,aes(x=class,y=distance))+geom_violin()+geom_boxplot(width=.12, fill=I('black'), notch=T, col='grey40')+ylab('SNPs distance at boundary region')+scale_x_discrete(labels=c("Misjoined haplotigs","Random error haplotigs"))+xlab(NULL)+coord_cartesian(ylim = c(0, 60000)) 
dev.off()

# size and location
par(mfrow=c(1,1))
#plot(random_error$n_end_position,random_error$distance,pch='.')
#plot(random_error$n_end_position,random_error$norm_distance,pch='.')
print('difference of median test for last snp or random error')
wilcox.test(distance~edge,data=random_error)
pdf(paste('box_plot_boundary_size_random_error_edge_effect','.pdf',sep=""))
ggplot(random_error,aes(x=edge,y=distance))+geom_violin()+geom_boxplot(width=.12, fill=I('black'), notch=T, col='grey40')+ylab('SNPs distance at boundary region')+scale_x_discrete(labels=c("Error at the first/last SNPs","Error at other SNPs"))+xlab(NULL)+coord_cartesian(ylim = c(0, 60000)) 
dev.off()
