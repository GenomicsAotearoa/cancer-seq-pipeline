
module load FastQC
mkdir -p fastqc_out
fastqc -t 1 $1 -o fastqc_out/
