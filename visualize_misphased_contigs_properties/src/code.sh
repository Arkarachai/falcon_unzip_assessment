#!/bin/bash
# visualize_misphased_contigs_properties 1.0.0

set -e -x -o pipefail

main() {

    echo "Value of filtered_informative_snps: '$filtered_informative_snps'"
    echo "Value of raw_plotting_files: '$raw_plotting_files'"

    dx download "$filtered_informative_snps" -o filtered_informative_snps


    ###################<execution start>###################

    ## generate all precalculate files
    python visual_falcon_unzip_prep.py filtered_informative_snps
    
    ## check longest contig for each type
    Rscript longest_track.R

    ## generating scatter_n_error_vs_switch_error
    Rscript scatter_n_error_vs_switch_error.R
    
    ## generating box_plot_boundary_size, histogram_normalize_boundary_location, box_plot_boundary_size_random_error_edge_effect
    Rscript boundary_visualizer.R 


    ###################<execution end>###################

    if [ "$raw_plotting_files" == true ]; then
        intermediate_file=$(dx upload aggregate_parent_switch_mis_join --brief)
        dx-jobutil-add-output intermediate_file "$intermediate_file" --class=array:file
        intermediate_file=$(dx upload aggregate_parent_switch_random_error --brief)
        dx-jobutil-add-output intermediate_file "$intermediate_file" --class=array:file
        intermediate_file=$(dx upload haplotig_evaluation_table --brief)
        dx-jobutil-add-output intermediate_file "$intermediate_file" --class=array:file
        intermediate_file=$(dx upload boundary_size_vs_location --brief)
        dx-jobutil-add-output intermediate_file "$intermediate_file" --class=array:file
    fi
    ls

    scatter_n_error_vs_switch_error=$(dx upload scatter_n_error_vs_switch_error.pdf --brief)
    box_plot_boundary_size=$(dx upload box_plot_boundary_size.pdf --brief)
    histogram_normalize_boundary_location=$(dx upload histogram_normalize_boundary_location.pdf --brief)
    box_plot_boundary_size_random_error_edge_effect=$(dx upload box_plot_boundary_size_random_error_edge_effect.pdf --brief)

    dx-jobutil-add-output scatter_n_error_vs_switch_error "$scatter_n_error_vs_switch_error" --class=file
    dx-jobutil-add-output box_plot_boundary_size "$box_plot_boundary_size" --class=file
    dx-jobutil-add-output histogram_normalize_boundary_location "$histogram_normalize_boundary_location" --class=file
    dx-jobutil-add-output box_plot_boundary_size_random_error_edge_effect "$box_plot_boundary_size_random_error_edge_effect" --class=file
}
