#!/bin/bash
#SBATCH --job-name=fastqc_combined
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --partition=general
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=user@uchc.edu
#SBATCH -o fastqc_combined_%j.out
#SBATCH -e fastqc_combined_%j.err
source ~/.bashrc
cd ~/DATA/RNA-seq/<Project_name>/merged_cf/combined_data/

fastqc -o ~/ANALYSIS/RNA-seq/<Project_name>/merged_cf/combined_data ~/DATA/RNA-seq/<Project_name>/merged_cf/combined_data/*.fastq*
