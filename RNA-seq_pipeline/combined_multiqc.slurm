#!/bin/bash
#SBATCH --job-name=combined_multiqc
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --partition=general
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=user@uchc.edu
#SBATCH -o combined_multiqc_%j.out
#SBATCH -e combined_multiqc_%j.err
source ~/.bashrc
cd ~/ANALYSIS/RNA-seq/<Project_name>/merged_cf/combined_data

multiqc -f -i RNA-seq_combined_fastq_report -n RNA-seq_combined_fastq_multiqc .