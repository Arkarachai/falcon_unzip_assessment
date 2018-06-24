# generate_ucsc_track_data (DNAnexus Platform App)

## What does this applet do?
This applet generates the UCSC genome browser custom track file(s) of mis-phased haplotigs from FALCON-Unzip algorithm based on trio Illumina data. It is provided as a reference to the publication **How well we can phase the assembled diploid human genome?: Assessment of FALCON-Unzip phasing using human trio** by *Arkarachai Fungtammasan* and *Brett Hannigan*. 

## What data are required for this applet to run?
The inputs include:

1 Table of filter haplotype and informative SNPs in tab-separated values format. This is the output of applet `falcon_unzip_evaluation`.

2 contigs id(s) of interest separated by comma

## What does this applet output?
The output includes custom track file(s) which can be view by UCSC genome browser.

## How does this applet work?
This applet is the shell script that runs a custom python scripts to create custom track according to UCSC genome browser specification. The underline python script can be retrieved using `dx get generate_ucsc_track_data`. The shell script which controls the python script will be stored in `src/code.sh`, while the python script will be stored in `resources/home/dnanexus`. 