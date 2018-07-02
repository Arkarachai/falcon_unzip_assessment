#!/bin/bash
# mummer_nucmer_aligner 1.0.0

set -x -o pipefail
main() {

    echo "Value of ref_fastagz: '$ref_fastagz'"
    echo "Value of query_fastagz: '$query_fastagz'"

    # The following line(s) use the dx command-line tool to download your file
    # inputs to the local file system using variable names for the filenames. To
    # recover the original filenames, you can use the output of "dx describe
    # "$variable" --name".
    refname=$(dx describe "$ref_fastagz" --name)
    refname="${refname%.gz}"
    refname="${refname%.fa}"
    refname="${refname%.fasta}"
    dx download "$ref_fastagz" -o "${refname}.gz"
    gunzip "${refname}.gz" &
    queryname=$(dx describe "$query_fastagz" --name)
    queryname="${queryname%.gz}"
    queryname="${queryname%.fa}"
    queryname="${queryname%.fasta}"
    dx download "$query_fastagz" -o "${queryname}.gz"
    gunzip "${queryname}.gz" &
    wait -n

    ###################<execution start>###################

    path_to_mummer="/MUMmer3.23"
    mkdir -p out/nucmer_full_outputs
    mkdir -p out/nucmer_filtered_outputs
    mkdir -p out/nucmer_filtered_snps
    
    if [ -n "$prefix" ]; then
        prefix="$prefix"
    else
        prefix="${queryname}_${refname}"
    fi
    
    ${path_to_mummer}/nucmer --prefix="$prefix" $extra_cmd "$refname" "$queryname"
    ${path_to_mummer}/show-coords -rcl "$prefix".delta > "$prefix".coords
    ${path_to_mummer}/delta-filter -r -q "$prefix".delta > "$prefix"_reduce.delta  
    ${path_to_mummer}/show-tiling "$prefix".delta > "$prefix".tiling
    ${path_to_mummer}/show-snps -Clr "$prefix".delta > "$prefix".snps
    ${path_to_mummer}/show-tiling "$prefix"_reduce.delta > "$prefix"_reduce.tiling
    ${path_to_mummer}/show-snps -Clr -T "$prefix"_reduce.delta > "$prefix"_reduce.snps
    mv "$prefix"_reduce.snps out/nucmer_filtered_snps
    mv "$prefix"_reduce.delta  "$prefix"_reduce.tiling out/nucmer_filtered_outputs
    mv *.tiling *.delta *.coords *.snps out/nucmer_full_outputs
    
    ###################<execution end>###################

    dx-upload-all-outputs
}
