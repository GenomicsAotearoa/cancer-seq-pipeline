// Trial pipeline
// Needs cutdapt, samtools, bwa, picard, java, varscan2, fastqc, bedtools. 
// Annovar script table_annovar.pl is a customised version of that provided with annovar, needs to replace the default version.
//      i.e. you need to install annovar and then replace table_annovar.pl with the copy found with this pipeline

concurrency=1000
bwaIndex="/nesi/project/uoa00595/data/hg19/ucsc.hg19.fasta"
picardLoc = "/scale_wlg_nobackup/pan_migration/apps/easybuild/RHEL6.3/westmere/software/picard/2.1.0/picard.jar"
hg19db = "/nesi/project/uoa00571/data/annovar/db/"
annovarDB = "/nesi/project/uoa00571/data/annovar/db"

qc = {
    output.dir="../data/fastqc"
    //exec "module load FastQC"
    multi "fastqc -t 1 $input1 -o fastqc/",
        "fastqc -t 1 $input2 -o fastqc/",
        "fastqc -t 1 $input3 -o fastqc/",
        "fastqc -t 1 $input4 -o fastqc/"
}


trim = {
    output.dir="/nesi/nobackup/uoa02606/data/intermediateFiles"
    // this assumes that the input is illumina, as we're not autodetecting across a variety of possible adapters,
    // which would be possible with trim_glaore, but that requires changes to trim_galores output methods. 
    // cores=0 detects how many cores are available and uses them.
    exec "cutadapt -a AGATCGGAAGAGC --pair-filter=any --cores=0 -q 16 --minimum-length=50 -o $output1 -p $output2 $input1 $input2","trim"
    //exec "cutadapt -a AGATCGGAAGAGC --pair-filter=any --cores=0 -q 16 --minimum-length=50 -o $output3 -p $output4 $input3 $input4","trim"
}



align = { 
    output.dir="/nesi/nobackup/uoa02606/data/intermediateFiles"
    // align the tumour sample and the 
    exec "bwa mem -t 12 -R '@RG\\tID:${SAMPLE}_02\\tSM:$SAMPLE\\tLB:${SAMPLE}_02\\tPL:ILLUMINA' $bwaIndex $input1 $input2 | samtools sort -@12 -O BAM -o $output.bam", "align"

}


markDuplicates = {
    output.dir="/nesi/nobackup/uoa02606/data/intermediateFiles"
    exec "java -jar  $picardLoc MarkDuplicates I=$input1.bam O=$output1.bam M=$output._dup_metrics.txt CREATE_INDEX=TRUE", "removeDuplicates"


}

recalibrate = {
    output.dir="/nesi/nobackup/uoa02606/data/intermediateFiles"
    exec """ ../bin/gatk/gatk --java-options "-Xmx16G" BaseRecalibrator -R $bwaIndex -I $input.bam -O $output.table --known-sites /scale_wlg_persistent/filesets/project/uoa02606/data/gatk/dbsnp_138.hg19.vcf --known-sites /scale_wlg_persistent/filesets/project/uoa02606/data/gatk/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf""", "recalibrate"
    exec """ ../bin/gatk/gatk --java-options "-Xmx16G" ApplyBQSR -R $bwaIndex -I $input.bam -bqsr $output.table -O $output.bam""", "recalibrate"

}

haplotypeCalling = {
    output.dir="/nesi/nobackup/uoa02606/data/intermediateFiles"
    exec """/scale_wlg_persistent/filesets/project/uoa02606/bin/gatk/gatk --java-options "-Xmx8G" HaplotypeCaller -R $bwaIndex -I $input.bam -O $output.vcf""" ,"haplotyping"
}


removeSuppplementary = {
    output.dir="/nesi/nobackup/uoa02606/data/intermediateFiles"
    exec "samtools view -@ 6 -b -F 2048 $input1.bam > $output1.bam", "removeSuppplementary"
    exec "samtools index $output1.bam", "removeSuppplementary"
    
}

pileUp = {
    output.dir="/nesi/nobackup/uoa02606/data/intermediateFiles"
    exec "samtools mpileup -P $threads -B -d 9000 -q 1 -f $bwaIndex  $input1.bam | java -Xmx16g -d64 -jar /nesi/project/uoa00571/bin/VarScan.v2.4.3.jar somatic --mpileup 1 --min-var-freq 0.1 --p-value 1.00 --somatic-p-value 1.00 --strand-filter 0 --tumor-purity 0.5 --output-vcf 1 --min-coverage-normal 10 --min-coverage-tumor 10 --output-snp $output1.snp.vcf  --output-indel $output1.indel.vcf","pileUp"

}


annotation = {
    output.dir="/nesi/nobackup/uoa02606/data/intermediateFiles"

    // Version 2 of our annotation
	exec """perl ./table_annovar.pl $input $annovarDB  -buildver hg19 -outfile $output.vcf -remove -protocol refGene,popfreq_max_20150413,dbnsfp30a,tfbsConsSites,targetScanS,genomicSuperDups,clinvar_20170130,cosmic82,avsnp147,rmsk,wgEncodeDacMapabilityConsensusExcludable,wgEncodeDukeMapabilityRegionsExcludable,gnomad_genome,gnomad_exome,exac03,intervar_20170202,revel,dgvMerged  -operation g,f,f,r,r,r,f,f,f,r,r,r,f,f,f,f,f,r  -nastring . -vcfinput""", "annotation"

}

count = {
  output.dir="/nesi/nobackup/uoa02606/data/intermediateFiles" 
  exec "cp $input1.vcf temp1.vcf"
  zipped = "$input1" + ".gz" 
  exec "bgzip -f $input1.vcf" 
  exec "tabix -f -p vcf $zipped"
  exec "cp temp1.vcf $input1.vcf"


  exec "cp $input2.vcf temp2.vcf"
  exec "bgzip -f $input2.vcf"
  zipped2 = "$input2" + ".gz"
  exec "tabix -f -p vcf $zipped2"
  exec "cp temp2.vcf $input2.vcf"


  //This really doesn't fall over gracefully when the filters are to strict.
  //Also, this is appallingly bad form. The only difference between the two calls are a bit flag and the presence/absence of some filters. 
  //it should be a single function called twice with different parameters.
  exec "Rscript src/filePairs.R $zipped $zipped2 $output"
  exec "Rscript src/filePairsSomatic.R $zipped $zipped2 $output2"

}

reportGeneration ={
    output.dir="../results"
   exec "Rscript src/generateReport.R $input1 $input2 $output1"
}

cleanup_bams = {
        cleanup "*.bam", "*.bam"
}


// currently runs qc, doeen't wait for confirmation on whether to continue the run, 
// trims, and does a 2nd lot of qc to the trimmed reads. 
// Once this is parrallelized, should be able to insert multiqc steps. 



run { 
   "%_*.fastq.gz" * [trim] + "%_*.trim" * [align] + "%_*" * [markDuplicates] + "%_*" * [recalibrate] + "%_" * [haplotypeCalling] +   "%_" *  [annotation]
}

