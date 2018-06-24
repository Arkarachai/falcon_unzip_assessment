#!/bin/bash
# falcon_unzip_evaluation 1.0.0

set -x -e -o pipefail
main() {

    ## download inputs
    echo "Value of trio_genotype_vcfgz: '$trio_genotype_vcfgz'"
    echo "Value of primary_contigs_nucmer: '$primary_contigs_nucmer'"
    echo "Value of haplotype_contigs_nucmer: '$haplotype_contigs_nucmer'"
    echo "Value of intermediate_results: '$intermediate_results'"

    dx download "$trio_genotype_vcfgz" -o trio_genotype_vcf.gz
    dx download "$primary_contigs_nucmer" -o primary_contigs_nucmer
    dx download "$haplotype_contigs_nucmer" -o haplotype_contigs_nucmer

    ###################<execution start>###################

    echo "### Filtering ###"
    ## pre-process primary and haplotype contig nucmer
    #1 get rid of header for primary and haplotype contig nucmer output
    tail -n +5 primary_contigs_nucmer > reduce_noheader_p.variants
    tail -n +5 haplotype_contigs_nucmer > reduce_noheader_h.variants
    number_haplotigs=$(cat reduce_noheader_h.variants | awk -F '\t' '{print $NF}' | sort | uniq | wc -l)
    echo "Haplotigs that map to reference genome with 1:1 relationship =" $number_haplotigs "haplotigs"  
    #2 get rid of parts of haplotype contigs are outside the corresponding primary contigs
    python mapping_range.py reduce_noheader_p.variants reduce_noheader_h.variants > reduce_noheader_h_pmatch.variants 
    number_haplotigs=$(cat reduce_noheader_h_pmatch.variants  | awk -F '\t' '{print $NF}' | sort | uniq | wc -l)
    echo "Haplotigs that mapped within primary contigs =" $number_haplotigs "haplotigs"  
    
    ## pre-process Illumina trio data
    #3 get all informative SNPs based on Illumina trio data
    zcat < trio_genotype_vcf.gz | python get_informative_snp.py > informative_snps
    
    ## Integrate Illumina with nucmer data
    #4 join haplotigs with informative variants
    python join_haplotype.py informative_snps reduce_noheader_h_pmatch.variants > informative_haplotype.snps
    #cat informative_haplotype.snps | cut -f 1 | sort | uniq -c
    #5 join primary contigs with informative variants
    python join_haplotype.py informative_snps reduce_noheader_p.variants > informative_primary.snps
    #cat informative_primary.snps | cut -f 1 | sort | uniq -c
    #6 interlace primary contig to each haplotig; now all SNPs on primary contigs will be in their corresponding haplotigs
    python interlace.py informative_primary.snps informative_haplotype.snps > interlace_p_h.snps
    #7 sort by chromosome and position
    cat interlace_p_h.snps | sort -k 12,12 -k 2n,2 > interlace_p_h_sorted.snps 
    number_snps=$(cat interlace_p_h_sorted.snps | wc -l )
    number_haplotigs=$(cut -f 13 interlace_p_h_sorted.snps | sort | uniq | wc -l)
    echo "Haplotigs that contain informative SNPs =" $number_snps "SNPs," $number_haplotigs "haplotigs"  
    #8 remove incorrect segregation between variants on haplotype contig and primary contig
    python removedup.py interlace_p_h_sorted.snps 2,4,12 > interlace_p_h_sorted_correctseg.snps
    number_snps=$(cat interlace_p_h_sorted_correctseg.snps | wc -l )
    number_haplotigs=$(cut -f 13 interlace_p_h_sorted_correctseg.snps | sort | uniq | wc -l)
    echo "Haplotigs with variants that report on either haplotig or primary contig =" $number_snps "SNPs," $number_haplotigs "haplotigs" 
    #9 remove haplotigs with less than 10 informative SNPs mapping.
    python filter_col_frequency.py interlace_p_h_sorted_correctseg.snps 13 10 > interlace_p_h_sorted_correctseg_min10.snps
    number_snps=$(cat interlace_p_h_sorted_correctseg_min10.snps | wc -l )
    number_haplotigs=$(cut -f 13 interlace_p_h_sorted_correctseg_min10.snps | sort | uniq | wc -l)
    echo "Haplotigs that have at least 10 informative SNPs =" $number_snps "SNPs," $number_haplotigs "haplotigs" 
    #10 remove haplotigs that span two chromosome
    python remove_haplotig_span_multiple_chr.py interlace_p_h_sorted_correctseg_min10.snps > interlace_p_h_sorted_correctseg_min10_1chr.snps
    number_snps=$(cat interlace_p_h_sorted_correctseg_min10_1chr.snps| wc -l )
    number_haplotigs=$(cut -f 13 interlace_p_h_sorted_correctseg_min10_1chr.snps | sort | uniq | wc -l)
    echo "Haplotigs that map to only one chromosome in reference genome =" $number_snps "SNPs," $number_haplotigs "haplotigs" 
    ##11 remove haplotigs that are overlapped other haplotigs
    #11 remove haplotigs that are interrupted by other haplotigs
    python remove_untillable_contig.py interlace_p_h_sorted_correctseg_min10_1chr.snps 13 > interlace_p_h_sorted_correctseg_tilled.snps
    #python remove_untillable_overlap_contig.py interlace_p_h_sorted_correctseg_min10_1chr.snps 13 > interlace_p_h_sorted_correctseg_tilled.snps
    number_snps=$(cat interlace_p_h_sorted_correctseg_tilled.snps | wc -l )
    number_haplotigs=$(cut -f 13 interlace_p_h_sorted_correctseg_tilled.snps | sort | uniq | wc -l)
    #echo "Haplotigs that do not overlap with other haplotigs =" $number_snps "SNPs," $number_haplotigs "haplotigs" 
    echo "Haplotigs that map continuously with reference genome =" $number_snps "SNPs," $number_haplotigs "haplotigs"
    #12 get rid of chrY and chrUn; note that this step should perform differently if other type of chromosome present in the data. In this case, there are only autosome, chrY, chrX and chrUn.
    cat  interlace_p_h_sorted_correctseg_tilled.snps | grep -v chrUn | grep -v chrY | grep -v chrM | grep -v chrEBV > interlace_p_h_sorted_correctseg_tilled_diploid1
    #13 keep only pseudo autosomal region for chrX
    python remove_x_specific.py interlace_p_h_sorted_correctseg_tilled_diploid1 > interlace_p_h_sorted_correctseg_tilled_diploid2
    number_snps=$(cat interlace_p_h_sorted_correctseg_tilled_diploid2 | wc -l )
    number_haplotigs=$(cut -f 13 interlace_p_h_sorted_correctseg_tilled_diploid2 | sort | uniq | wc -l)
    echo "Haplotigs that map to autosome and pseudoautosomal regions =" $number_snps "SNPs," $number_haplotigs "haplotigs" 
    echo ""

    ## Data analysis
    #14 generate haplotig_evaluation_table
    python error_rate_prep.py interlace_p_h_sorted_correctseg_tilled_diploid2 > haplotig_evaluation_table 
    echo "### Contig type ###"
    cat haplotig_evaluation_table | awk -F '\t' '{print $NF}' | tail -n +2 | sort | uniq -c
    echo ''
    echo "### Primvary vs haplotig contribution ###"
    cat interlace_p_h_sorted_correctseg_tilled_diploid2 | awk -F '\t' '{print $NF}' | sort | uniq -c
    echo ''
    echo "### Paternal vs maternal contribution ###"
    cat interlace_p_h_sorted_correctseg_tilled_diploid2 | awk -F '\t' '{print $1}' | sort | uniq -c
    
    ###################<execution end>###################

    ## upload outputs
    informative_snps=$(dx upload informative_snps --brief)
    filtered_informative_snps=$(dx upload interlace_p_h_sorted_correctseg_tilled_diploid2 -o filtered_informative_snps --brief)
    haplotig_evaluation_table=$(dx upload haplotig_evaluation_table --brief)

    dx-jobutil-add-output informative_snps "$informative_snps" --class=file
    dx-jobutil-add-output filtered_informative_snps "$filtered_informative_snps" --class=file
    dx-jobutil-add-output haplotig_evaluation_table "$haplotig_evaluation_table" --class=file
    
    if "$intermediate_results" == true; then
        intermediate_results=$(dx upload reduce_noheader_p.variants --brief)
        dx-jobutil-add-output intermediate_results "$intermediate_results" --class=array:file
        intermediate_results=$(dx upload reduce_noheader_h.variants --brief)
        dx-jobutil-add-output intermediate_results "$intermediate_results" --class=array:file
        intermediate_results=$(dx upload reduce_noheader_h_pmatch.variants  --brief)
        dx-jobutil-add-output intermediate_results "$intermediate_results" --class=array:file
        intermediate_results=$(dx upload informative_haplotype.snps --brief)
        dx-jobutil-add-output intermediate_results "$intermediate_results" --class=array:file
        intermediate_results=$(dx upload informative_primary.snps --brief)
        dx-jobutil-add-output intermediate_results "$intermediate_results" --class=array:file
        intermediate_results=$(dx upload interlace_p_h.snps --brief)
        dx-jobutil-add-output intermediate_results "$intermediate_results" --class=array:file
        intermediate_results=$(dx upload interlace_p_h_sorted_correctseg.snps --brief)
        dx-jobutil-add-output intermediate_results "$intermediate_results" --class=array:file
        intermediate_results=$(dx upload interlace_p_h_sorted_correctseg_min10.snps --brief)
        dx-jobutil-add-output intermediate_results "$intermediate_results" --class=array:file
        intermediate_results=$(dx upload interlace_p_h_sorted_correctseg_min10_1chr.snps --brief)
        dx-jobutil-add-output intermediate_results "$intermediate_results" --class=array:file
        intermediate_results=$(dx upload interlace_p_h_sorted_correctseg_tilled.snps --brief)
        dx-jobutil-add-output intermediate_results "$intermediate_results" --class=array:file
        intermediate_results=$(dx upload interlace_p_h_sorted_correctseg_tilled_diploid2 --brief)
        dx-jobutil-add-output intermediate_results "$intermediate_results" --class=array:file
    fi
}
