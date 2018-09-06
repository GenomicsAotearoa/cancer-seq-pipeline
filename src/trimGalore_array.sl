#!/bin/bash
#SBATCH -A uoa00571
#SBATCH --job-name       TrimGalore_array
#SBATCH --cpus-per-task  1
#SBATCH --time           6:00:00
#SBATCH --mem            8G
#SBATCH -o trimGalore_%a.out # Standard output
#SBATCH -e trimGalore_%a.err # Standard error

module load TrimGalore/0.4.2-foss-2015a

R1FILES=($(ls -1 raw/*_1.fastq.gz))

R1=${R1FILES[$SLURM_ARRAY_TASK_ID]}
R2=${R1%_1.fastq.gz}_2.fastq.gz

mkdir -p trimmed_reads # out directory

trim_galore -q 30 --length 50 --paired -o trimmed_reads $R1 $R2
