# falcon_unzip_evaluation (DNAnexus Platform App)

## What does this applet do?
This applet evaluates of the FALCON-Unzip algorithm using trio Illumina data. It is provided as a reference to the publication **How well we can phase the assembled diploid human genome?: Assessment of FALCON-Unzip phasing using human trio** by *Arkarachai Fungtammasan* and *Brett Hannigan*. It aims for reproducibility. The users are free to make use of all the underline scripts for their own researches related/unrelated to this study. However, with a different experimental set including changing the version of reference genome or using a daughter instead of a son, some part of the scripts should be changed to reflect different experimental designs.

## What data are required for this applet to run?
The inputs include:

1 Multi-sample VCF of the trio. The order of trio go by son-father-mother

2 Variants calling of assembled son primary contigs from MUMmur program

3 Variants calling of assembled son haplotype contigs from MUMmur program

## What does this applet output?
The outputs include:

1 Informative SNPs, SNPs that can be determined parental origin, from trio data in tab-separated values format. Each row is one SNP and it represent:
`allele_inherit_from_father` `allele_inherit_from_mother` `CHROM` `POS` `ID` `REF` `ALT`     `QUAL` `FILTER` `INFO` `FORMAT` `son` `father` `mother` `trio_allele_code`

2 Table of filter haplotype and informative SNPs in tab-separated values format. Each row is one SNP and it represents:
`parent_of_origin_for_SNP_either_dad_or_mom` 
`position_of_SNP_in_reference`
`character_at_this_position_in_reference`
`character_at_this_position_in_query`
`position_of_SNP_in_query`
`distance_from_this_SNP_to_nearest_mismatch`
`distance_from_this_SNP_to_nearest_sequence_end`
`length_of_the_reference_sequence`
`length_of_the_query_sequence`
`sequence_direction`
`the_reference_FastA_IDs`
`query_FastA_IDs`
`contig_of_origin_for_SNPs_either_primary_or_haplotig`

The second to the thirteen column are the outputs from show-snps function of mummer
http://mummer.sourceforge.net/manual/

3 Haplotype evaluation table. Each row is one haplotype it represents:
`CHROM` `contig_name` `number_informative_snps` `error_rate` `snp_span_length` `parent_switch_n` `contig_origin_switch_n` `parent_pattern` `contig_origin_pattern`
`pattern_switch_position` `contig_origin_switch_position` `parent_contig_subset`
`haplotig_daddy_rate` `error_string` `error_count` `count_i` `count_ii` `count_iii` `count_iiii` `count_i>5` `phasing_type`

## How does this applet work?
This applet is the shell script that ties together several custom python scripts to filter and evaluate phasing accuracy of FALCON-Unzip algorithm. The underline python scripts can be retrieved using `dx get falcon_unzip_evaluation`. The shell script which orchestrates the analysis will be stored in `src/code.sh`, while all the python scripts will be stored in `resources/home/dnanexus`. 

