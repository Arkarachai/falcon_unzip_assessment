## **How well we can phase the assembled diploid human genome?: Assessment of FALCON-Unzip phasing using human trio**

The purpose of this GitHub repo is to distribute the analysis code of the publication **How well we can phase the assembled diploid human genome?: Assessment of FALCON-Unzip phasing using human trio** by *Arkarachai Fungtammasan* and *Brett Hannigan*. It aims for reproducibility. The users are free to make use of all the underline scripts for their own researches related/unrelated to this study. However, with a different experimental set including changing the version of reference genome or using a daughter instead of a son, some part of the scripts should be changed to reflect different experimental designs.

These applets have already been built to the public project-FGvfpX801xQ9Q4GZ6FfPZGxy. However, the user may also git clone this GitHub and build them yourself.

Although the code is organized to operate on DNAnexus, the user can also run them on their local computer. 

1 Please install required binary to your local computer. Each tool has the following dependency.
1. mummer_nucmer_aligner: MUMmer3(at root), numpy
2. falcon_unzip_evaluation: None
3. visualize_misphased_contigs_properties: ggplot2
4. generate_ucsc_track_data: python-pandas

The dependency information is also stored in the asset/Makefile for each tool folder, except MUMmer3 which has to be download from the website directly.

2 Copy the script from resources/home/dnanexus/* to working directory.

3 Run src/code.sh between `<execution start>` and `<execution end>`.
