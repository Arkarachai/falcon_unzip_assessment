dataset=read.table('haplotig_evaluation_table',sep='\t',header=T)

no_error=dataset[dataset$phasing_type=='no_error',]
all_error=dataset[dataset$phasing_type!='no_error',]
lowset=dataset[dataset$phasing_type=='random_error',]
highset=dataset[dataset$phasing_type=='mis_join',]
#mediumset=data.frame(rbind(dataset[dataset$phasing_type=='inbetween',],dataset[dataset$phasing_type=='mix',]))
#mixset=dataset[dataset$phasing_type=='mix',]
#inbetweenset=dataset[dataset$phasing_type=='inbetween',]

print('no_error_stat')
paste('total_area_no_error',sum(no_error$snp_span_length))
paste('number of contig',nrow(no_error))
paste('longest_no_error', no_error[no_error$snp_span_length==max(no_error$snp_span_length),]$contig,'contain',
      max(no_error$snp_span_length),'bp',no_error[no_error$snp_span_length==max(no_error$snp_span_length),]$N,'snps')

print('error_stat')
paste('total_area_error',sum(all_error$snp_span_length))
paste('number of contig',nrow(all_error))
paste('max_error=',max(all_error$error_rate))
paste('min_error=',min(all_error$error_rate))
paste('single_misphased_contig', nrow(all_error[all_error$error_count==1,]),'contigs; covering',
      sum(all_error[all_error$error_count==1,]$snp_span_length))
paste('single_misphased_contig', nrow(all_error[all_error$error_rate<0.10,]),'contigs; covering',
      sum(all_error[all_error$error_rate<0.10,]$snp_span_length))


paste('longest_all_error', all_error[all_error$snp_span_length==max(all_error$snp_span_length),]$contig,'contain',
      max(all_error$snp_span_length),'bp',all_error[all_error$snp_span_length==max(all_error$snp_span_length),]$N,'snps')

paste('longest_misjoin', highset[highset$snp_span_length==max(highset$snp_span_length),]$contig,'contain',
      max(highset$snp_span_length),'bp',highset[highset$snp_span_length==max(highset$snp_span_length),]$N,'snps')

paste('longest_random_error', lowset[lowset$snp_span_length==max(lowset$snp_span_length),]$contig,'contain',
      max(lowset$snp_span_length),'bp',lowset[lowset$snp_span_length==max(lowset$snp_span_length),]$N,'snps')



print('total_area_random_error')
sum(lowset$snp_span_length)
print('total_area_mis_join')
sum(highset$snp_span_length)


paste('highest_snp_density_no_error',no_error[no_error$N==max(no_error$N),]$contig)
paste('highest_snp_density_random_error',lowset[lowset$N==max(lowset$N),]$contig)
paste('highest_snp_density_mis_join',highset[highset$N==max(highset$N),]$contig)

paste('highest_error_count_random_error',lowset[lowset$error_count==max(lowset$error_count),]$contig)
paste('highest_error_count_mis_join',highset[highset$error_count==max(highset$error_count),]$contig)
paste('highest_error_count_all',all_error[all_error$error_count==max(all_error$error_count),]$contig)

paste('highest_error_rate_random_error',lowset[lowset$error_rate==max(lowset$error_rate),]$contig)
paste('highest_error_rate_mis_join',highset[highset$error_rate==max(highset$error_rate),]$contig)
paste('highest_error_rate_all',all_error[all_error$error_rate==max(all_error$error_rate),]$contig)
