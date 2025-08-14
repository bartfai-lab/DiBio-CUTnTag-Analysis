DiBioCUT&Tag Analysis

This project aims to provide the exact code used for the analysis of our next generation sequencing data in our publication "CUT&Tag and DiBioCUT&Tag enable investigation of the AT-rich epigenome of Plasmodium falciparum from low input samples".

If you want to recreate our datasets or apply our scripts to your own data this repository can provide you with the recources you need. 
CT_mapping_bowtie2v4.sh is a script to generate bigWig / bedgraph files for visualisation of CUT&Tag datasets on e.g UCSC genome browser from paired end sequencing data, CT_mapping_bowtie2v4.singleend.sh is the same script for single end data)
CT_log2_ratio_binsv2.sh is a script to normalise CUT&Tag tracks to a control using a log2 ratio normalisation. Due to the sparce nature of CUT&Tag tracks, data is first binned.
All raw and processed sequencing files from our study have been submitted to Gene Expression Omnibus (GEO) under the reference number GSE270104.

To get started with our CUT&Tag mapping script (CT_mapping_bowtie2v4.sh):
1. set up a new conda environment with the provided .yml file CTanalysis.yml.
2. index your reference genome with bowtie2
3. Generate a working directory containing a /fastq folder with the fastq files you want to analyse.
4. Run the script: bash /path/to/CT_mapping_bowtie2v4.sh -n name_of_paired_fastq_without_R*.fastq.gz -g /path/to/genome -p /path/to/project_directory -q desired_mapping_quality_cutoff
important: -p /path/to/project_directory needs to contain the /fastq folder for this script to work


To get started with our log2 ratio track generation script (CT_log2_ratio_bins.sh):
1. make sure you have deeptools installed and already generated bigWig files for the datasets you want to compare (e.g by running CT_mapping_bowtie2v4.sh)
2. run the script: bash /path/to/CT_log2_calc_binv2.sh -c /path/to/cutntag.bw  -b /path/to/background.bw -s desired_binsize -o /path/to/output -n outputfile_name


To get started with our bedgraph averaging script (BDG_averagesv2.sh):
1. make sure you have bedtools installed and already generated bedgraph files for the datasets you want to combine/average (e.g by running CT_mapping_bowtie2v4.sh)
2. run the script: bash /path/to/BDG_averagesv1.sh -g /path/to/chrNameLength.txt -p /path/to/project_directory -a /path/to/file1.bedgraph -b /path/to/file2.bedgraph -c /path/to/file3.bedgraph -n output_name 


If you need help with analysis of CUT&Tag data, we would like to refer you to the very comprehensive tutorial from Zheng Y. et al https://yezhengstat.github.io/CUTTag_tutorial/ which has heavily influenced our scripts.

This project is part of a publication and will therefore not be regularly updated.
