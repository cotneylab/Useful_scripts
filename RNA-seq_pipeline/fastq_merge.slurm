#!/bin/bash
#SBATCH --job-name=fastq_merge
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --partition=general
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=user@uchc.edu
#SBATCH -o fastq_merge_%j.out
#SBATCH -e fastq_merge_%j.err

source ~/.bashrc

# Merges Fastq files from each lane into one file for each read
# name given to the merged file can include the sample name, an identifier and the sample number (S1, S2, etc) assigned by illumina

# CS13-2_S14_L001_R1_001.fastq.gz
# CS13-2_S14_L002_R1_001.fastq.gz
# CS13-2_S14_L003_R1_001.fastq.gz
# CS13-2_S14_L004_R1_001.fastq.gz
# merge to AWIL_017_CS13-2_S14_R1.fastq
# where AWIL_017 is the sample name, CS13-2 is a sample identifier and S14 is the sample number

cd ~/DATA/RNA-seq/<Project_name>/
mkdir merged_cf

# Merge the Fastq files for each lane of R1 together, R2 together 
cd ~/DATA/RNA-seq/<Project_name>/<Sample_name>/Files
zcat *R1*fastq.gz > ~/DATA/RNA-seq/<Project_name>/merged_cf/<Sample_name>_<identifier>_<Sample_number>_R1.fastq
cd ~/DATA/RNA-seq/<Project_name>/<Sample_name>/Files
zcat *R2*fastq.gz > ~/DATA/RNA-seq/<Project_name>/merged_cf/<Sample_name>_<identifier>_<Sample_number>_R2.fastq
cd ~/DATA/RNA-seq/<Project_name>/merged_cf
gzip *.fastq
