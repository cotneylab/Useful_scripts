#!/bin/bash
#SBATCH --job-name=merged_fastqc
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --partition=general
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=user@uchc.edu
#SBATCH -o merged_fastqc_%j.out
#SBATCH -e merged_fastqc_%j.err

source ~/.bashrc

# Perform Fastqc on merged fastq files
cd ~/DATA/RNA-seq/<Project_name>/merged_cf
fastqc -o ~/ANALYSIS/RNA-seq/<Project_name>/ ~/DATA/RNA-seq/<Project_name>/merged_cf/*.fastq*
