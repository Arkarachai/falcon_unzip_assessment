#!/bin/bash
# generate_ucsc_track_data 1.0.0

set -e -x -o pipefail
main() {

    echo "Value of filtered_informative_snps: '$filtered_informative_snps'"
    echo "Value of contig_id: '$contig_id'"

    dx download "$filtered_informative_snps" -o filtered_informative_snps

    ###################<execution start>###################

    python ucsc_track_generator.py filtered_informative_snps "$contig_id"

    ###################<execution end>###################

    for i in $(ls *.ucsc.track); do
        custom_track_files=$(dx upload $i --brief)
        dx-jobutil-add-output custom_track_files "${custom_track_files}" --class=array:file
    done

}
