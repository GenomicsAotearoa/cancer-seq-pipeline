## Analysis pipeline

Changed repository name - this pipeline was initially being put together as part of the prosper pipeline, it is now being developed as a tool across multiple projects. 
Version 1 of bpipe pipeline for analysis of genomic data. 

# Requirements
Bpipe is installed (link to latest version available at docs.bpipe.org)
Tools that bpipe uses:
* FastQC
* cutadapt
* BWA
* SAMtools
* picard
* VCFtools
* BEDTools
* Perl
* R

Currently using a module system where tools need to be loaded prior to the execution of each step. 

# Running

Assuming there is a data directory at the same level as your src directory, then execute by moving to the data directory, and running:

bpipe ../src/pipeline tumour_1.fastq.gz tumour_2.fastq.gz normal_1.fastq.gz normal_2.fastq.gz

# Current Steps


# Todo
The report generation step requires a file in the data directory called metadata.txt, containing lists of genes. Currently, a single list is selected and filtered for during report generation.
Awaiting confirmation that the current list can go public. 

