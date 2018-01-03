#!/bin/bash
#SBATCH --job-name=merged_multiqc
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --partition=general
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=user@uchc.edu
#SBATCH -o merged_multiqc_%j.out
#SBATCH -e merged_multiqc_%j.err

source ~/.bashrc

# Perform multiqc on merged fastq files
cd ~/ANALYSIS/RNA-seq/<Project_name>
multiqc -f -i RNA-seq_merged_file_report -n RNA-seq_merged_multiqc .
