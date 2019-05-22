## Analysis pipeline

Changed repository name - this pipeline was initially being put together as part of the prosper pipeline, it is now being developed as a tool across multiple projects. 
Version 1 of bpipe pipeline for analysis of genomic data. 

# Requirements
Bpipe is installed (link to latest version available at docs.bpipe.org)
Tools that bpipe uses:
* [bpipe](https://github.com/ssadedin/bpipe) - also requires Java
* FastQC
* cutadapt
* BWA
* SAMtools
* picard
* VCFtools
* VarScan2
* BEDTools
* Perl
* R

Currently using a module system where tools need to be loaded prior to the execution of each step. 

# Running

Assuming there is a data directory at the same level as your src directory, then execute by moving to the data directory, and running:

```
bpipe ../src/pipeline tumour_1.fastq.gz tumour_2.fastq.gz normal_1.fastq.gz normal_2.fastq.gz
```

# Current Steps


# Todo/notes
The report generation step requires a file in the data directory called metadata.txt, containing lists of genes. Currently, a single list is selected and filtered for during report generation.
Awaiting confirmation that the current list can go public. 

Command to run: 

```
bpipe src/pipeline data/normal_1.fq.gz data/normal_2.fq.gz data/tumour_1.fq.gz data/tumour_2.fq.gz
```

Currently picks up the sample name from the name of the first file input. It picks up everything prior to the first underscore and us
es that. i.e. in the example above, the final report would list the filename as "normal"

Location of hg19 genome is currently hardcoded in another projects directory, i.e. the TestData directory. Not a problem, just something to be aware of/not the best form.

Varscan step is not currently operational. Needs a local version of varscan, which I've probably been told the location of but can't currently find/recall.

Currently running only serially. Not a problem for the otago box.

Need a place to centrally locate filter parameters. i.e. there should be a file where it can be checked easily, currently filter levels are buried in R code of filePair.R and generateReport.R

`filePairs.R` and `filePairsSomatic.R` really need to be condensed to a single step. The only difference between them is a flag stating somatic/germline and some filter settings. Also results in the VCF file's being read in twice, which is inelegant.

add the cleanup module to the pipeline to get rid of all the intermediate bams/save a ton of space.

`data/dataColumns.csv` lists the fields colallted by the filePairs scripts.

`data/reportColumns.csv` contains the fields that are automatically put into the final report. Though there are currently a couple in there that aren't, as I need to confirm with Cris that I'm putting the correct thing in there.

filePairs files should probably be renamed to something more reflective of function.

some of the files being generated by the R scripts are not being located properly, resulting in the untidy top level directory. This
needs to be fixed. i.e. the final resulting html file is not going to the results directory. 

Add some specific instructions on how to include the insert size selection AWK scripts

Add mutational signatures


