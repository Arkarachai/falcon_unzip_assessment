# mummer_nucmer_aligner (DNAnexus Platform App)

## What does this applet do?
This applet performs nucmer alignment between query sequences with the reference sequence. It is provided as a reference to the publication **How well we can phase the assembled diploid human genome?: Assessment of FALCON-Unzip phasing using human trio** by *Arkarachai Fungtammasan* and *Brett Hannigan*. 

## What data are required for this applet to run?
The inputs include:

1 Query FASTA sequence. In this study, we use primary `cns_p_ctg.fasta.gz` and haplotype contigs `cns_h_ctg.fasta.gz`

2 Reference FASTA sequence. In this study, we used `GRCh38.no_alt_analysis_set.fa.gz`


## What does this applet output?
The outputs include:

1 The native nucmer results using function `nucmer`, `show-coord`, `show-tiling`, and `show-snps`

2 The `show-snps`, `show-coord` and `show-tiling` outputs from scope search of `nucmer` alignment using `delta-filter -r -q`. 


## How does this applet work?
This applet is the shell script that orchestrates MUMmer3 suite. The underline MUMmer3 can be downloaded from mummer.sourceforge.net. The shell script which orchestrates the MUMmer3 will be stored in `src/code.sh` after the user runs `dx get mummer_nucmer_aligner`.