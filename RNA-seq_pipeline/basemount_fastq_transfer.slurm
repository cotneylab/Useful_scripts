#!/bin/bash
#SBATCH --job-name=Basemount_fastq_transfer
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --partition=general
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=user@uchc.edu
#SBATCH -o Basemount_fastq_transfer_%j.out
#SBATCH -e Basemount_fastq_transfer_%j.err

source ~/.bashrc
module load basemount

cd ~
mkdir -p /tmp/user/basespace
basemount /tmp/user/basespace

cd /tmp/user/basespace

cp -r Projects/<Project_name>/Samples/<Sample_name> ~/DATA/RNA-seq/<Project_name>/

cd ~
basemount --unmount /tmp/user/basespace
rm -r tmp/user/
