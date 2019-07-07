#!/bin/bash
# falcon_unzip_evaluation 1.0.0

set -x -e -o pipefail
main() {

    ## download inputs
    echo "Value of trio_genotype_vcfgz: '$trio_genotype_vcfgz'"
    echo "Value of primary_contigs_nucmer: '$primary_contigs_nucmer'"
    echo "Value of intermediate_results: '$intermediate_results'"

    dx download "$trio_genotype_vcfgz" -o trio_genotype_vcf.gz
    dx download "$primary_contigs_nucmer" -o primary_contigs_nucmer

    ###################<execution start>###################

    echo "### Filtering ###"
    ## pre-process primary and haplotype contig nucmer
    #1 get rid of header for primary and haplotype contig nucmer output
    tail -n +5 primary_contigs_nucmer > reduce_noheader_p.variants
    number_haplotigs=$(cat reduce_noheader_p.variants | awk -F '\t' '{print $NF}' | sort | uniq | wc -l)
    echo "Haplotigs that map to reference genome with 1:1 relationship =" $number_haplotigs "haplotigs"  
       
    ## pre-process Illumina trio data
    #2 get all informative SNPs based on Illumina trio data
    zcat < trio_genotype_vcf.gz | python get_informative_snp.py > informative_snps
    
    ## Integrate Illumina with nucmer data
    #cat informative_haplotype.snps | cut -f 1 | sort | uniq -c
    #3 join primary contigs with informative variants
    python join_haplotype.py informative_snps reduce_noheader_p.variants > all_informative_primary.snps
    cat all_informative_primary.snps | grep -v ^problematic > informative_primary.snps
    echo "snps origin"
    cat informative_primary.snps | cut -f 1 | sort | uniq -c
    #4 sort by chromosome and position
    cat informative_primary.snps | sort -k 12,12 -k 2n,2 > interlace_p_h_sorted_sort.snps 
    number_snps=$(cat interlace_p_h_sorted_sort.snps | wc -l )
    number_haplotigs=$(cut -f 13 interlace_p_h_sorted_sort.snps | sort | uniq | wc -l)
    echo "Haplotigs that contain informative SNPs =" $number_snps "SNPs," $number_haplotigs "haplotigs"  
    echo "Haplotigs with variants that report on either haplotig or primary contig =" $number_snps "SNPs," $number_haplotigs "haplotigs" 
    #5 remove haplotigs with less than 10 informative SNPs mapping.
    python filter_col_frequency.py interlace_p_h_sorted_sort.snps 13 10 > interlace_p_h_sorted_correctseg_min10.snps
    number_snps=$(cat interlace_p_h_sorted_correctseg_min10.snps | wc -l )
    number_haplotigs=$(cut -f 13 interlace_p_h_sorted_correctseg_min10.snps | sort | uniq | wc -l)
    echo "Haplotigs that have at least 10 informative SNPs =" $number_snps "SNPs," $number_haplotigs "haplotigs" 
    #6 remove haplotigs that span two chromosome
    python remove_haplotig_span_multiple_chr.py interlace_p_h_sorted_correctseg_min10.snps > interlace_p_h_sorted_correctseg_min10_1chr.snps
    number_snps=$(cat interlace_p_h_sorted_correctseg_min10_1chr.snps| wc -l )
    number_haplotigs=$(cut -f 13 interlace_p_h_sorted_correctseg_min10_1chr.snps | sort | uniq | wc -l)
    echo "Haplotigs that map to only one chromosome in reference genome =" $number_snps "SNPs," $number_haplotigs "haplotigs" 
    ##7 remove haplotigs that are overlapped other haplotigs
    python remove_untillable_contig.py interlace_p_h_sorted_correctseg_min10_1chr.snps 13 > interlace_p_h_sorted_correctseg_tilled.snps
    #python remove_untillable_overlap_contig.py interlace_p_h_sorted_correctseg_min10_1chr.snps 13 > interlace_p_h_sorted_correctseg_tilled.snps
    number_snps=$(cat interlace_p_h_sorted_correctseg_tilled.snps | wc -l )
    number_haplotigs=$(cut -f 13 interlace_p_h_sorted_correctseg_tilled.snps | sort | uniq | wc -l)
    #echo "Haplotigs that do not overlap with other haplotigs =" $number_snps "SNPs," $number_haplotigs "haplotigs" 
    echo "Haplotigs that map continuously with reference genome =" $number_snps "SNPs," $number_haplotigs "haplotigs"
    #8 get rid of chrY and chrUn; note that this step should perform differently if other type of chromosome present in the data. In this case, there are only autosome, chrY, chrX and chrUn.
    cat  interlace_p_h_sorted_correctseg_tilled.snps | awk '$12!=Un' | awk '$12!=Y' | awk '$12!=M' | awk '$12!=EBV' > interlace_p_h_sorted_correctseg_tilled_diploid1
    #9 keep only pseudo autosomal region for chrX
    python remove_x_specific.py interlace_p_h_sorted_correctseg_tilled_diploid1 > interlace_p_h_sorted_correctseg_tilled_diploid2
    number_snps=$(cat interlace_p_h_sorted_correctseg_tilled_diploid2 | wc -l )
    number_haplotigs=$(cut -f 13 interlace_p_h_sorted_correctseg_tilled_diploid2 | sort | uniq | wc -l)
    echo "Haplotigs that map to autosome and pseudoautosomal regions =" $number_snps "SNPs," $number_haplotigs "haplotigs" 
    echo ""

    ## Data analysis
    #10 generate haplotig_evaluation_table
    python error_rate_prep.py interlace_p_h_sorted_correctseg_tilled_diploid2 $parent > haplotig_evaluation_table 
    echo "### accuracy type ###"
    cat haplotig_evaluation_table | awk -F '\t' '{print $NF}' | tail -n +2 | sort | uniq -c
    echo "### spanning snps ###"
    cat haplotig_evaluation_table | tail -n +2 | awk '{sum+=$5}; END {printf "%.3f\n", sum}'

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
    
    ls
    ## check longest contig for each type
    Rscript longest_track.R

    ## generating scatter_n_error_vs_switch_error
    #Rscript scatter_n_error_vs_switch_error.R
    Rscript scatter_n_error_vs_switch_error_all_type.R
    
    ## generating box_plot_boundary_size, histogram_normalize_boundary_location, box_plot_boundary_size_random_error_edge_effect
    #Rscript boundary_visualizer.R 

    ## generating switchplot for mis_join and random error
    #Rscript switch_plot.R aggregate_parent_switch_mis_join
    #Rscript switch_plot.R aggregate_parent_switch_random_error


    ###################<execution end>###################

    if [ "$intermediate_results" == true ]; then
        intermediate_results=$(dx upload reduce_noheader_p.variants --brief)
        dx-jobutil-add-output intermediate_results "$intermediate_results" --class=array:file
        intermediate_results=$(dx upload informative_primary.snps --brief)
        dx-jobutil-add-output intermediate_results "$intermediate_results" --class=array:file
        intermediate_results=$(dx upload interlace_p_h_sorted_correctseg_min10.snps --brief)
        dx-jobutil-add-output intermediate_results "$intermediate_results" --class=array:file
        intermediate_results=$(dx upload interlace_p_h_sorted_correctseg_min10_1chr.snps --brief)
        dx-jobutil-add-output intermediate_results "$intermediate_results" --class=array:file
        intermediate_results=$(dx upload interlace_p_h_sorted_correctseg_tilled.snps --brief)
        dx-jobutil-add-output intermediate_results "$intermediate_results" --class=array:file
        intermediate_results=$(dx upload interlace_p_h_sorted_correctseg_tilled_diploid2 --brief)
        dx-jobutil-add-output intermediate_results "$intermediate_results" --class=array:file
        intermediate_results=$(dx upload aggregate_parent_switch_mis_join --brief)
        dx-jobutil-add-output intermediate_results "$intermediate_results" --class=array:file
        intermediate_results=$(dx upload aggregate_parent_switch_random_error --brief)
        dx-jobutil-add-output intermediate_results "$intermediate_results" --class=array:file
        intermediate_results=$(dx upload boundary_size_vs_location --brief)
        dx-jobutil-add-output intermediate_results "$intermediate_results" --class=array:file
 #       intermediate_results=$(dx upload aggregate_parent_switch_mis_join.pdf --brief)
 #       dx-jobutil-add-output intermediate_results "$intermediate_results" --class=array:file
 #       intermediate_results=$(dx upload aggregate_parent_switch_random_error.pdf --brief)
 #       dx-jobutil-add-output intermediate_results "$intermediate_results" --class=array:file     
    fi

#    scatter_n_error_vs_switch_error=$(dx upload scatter_n_error_vs_switch_error.pdf --brief)
#    scatter_n_error_vs_switch_error_all_type=$(dx upload scatter_n_error_vs_switch_error_all_type.pdf --brief)
#    scatter_error_density_vs_switch_error_all_type=$(dx upload scatter_error_density_vs_switch_error_all_type.pdf --brief)
    scatter_error_ratio_vs_switch_error_all_type=$(dx upload scatter_error_ratio_vs_switch_error_all_type.pdf --brief)
 #   box_plot_boundary_size=$(dx upload box_plot_boundary_size.pdf --brief)
 #   histogram_normalize_boundary_location=$(dx upload histogram_normalize_boundary_location.pdf --brief)
 #   box_plot_boundary_size_random_error_edge_effect=$(dx upload box_plot_boundary_size_random_error_edge_effect.pdf --brief)

 #   dx-jobutil-add-output scatter_n_error_vs_switch_error_all_type "$scatter_n_error_vs_switch_error_all_type" --class=file
 #   dx-jobutil-add-output scatter_n_error_vs_switch_error "$scatter_n_error_vs_switch_error" --class=file
 #   dx-jobutil-add-output scatter_error_density_vs_switch_error_all_type "$scatter_error_density_vs_switch_error_all_type" --class=file
    dx-jobutil-add-output scatter_error_ratio_vs_switch_error_all_type "$scatter_error_ratio_vs_switch_error_all_type" --class=file
#    dx-jobutil-add-output box_plot_boundary_size "$box_plot_boundary_size" --class=file
#    dx-jobutil-add-output histogram_normalize_boundary_location "$histogram_normalize_boundary_location" --class=file
#    dx-jobutil-add-output box_plot_boundary_size_random_error_edge_effect "$box_plot_boundary_size_random_error_edge_effect" --class=file

}
