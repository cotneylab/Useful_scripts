echo "#!/bin/bash
#SBATCH --job-name=zcat_reps
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --partition=general
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=user@uchc.edu
#SBATCH -o zcat_reps_%j.out
#SBATCH -e zcat_reps_%j.err
source ~/.bashrc
cd ~/DATA/RNA-seq/<Project_name>/merged_cf" > zcat_reps.slurm

cat ~/DATA/RNA-seq/<Project_name>/merged_cf/combined_data/sample_list.txt | awk '{
print "zcat "$1"*R1*fastq > ~/DATA/RNA-seq/<Project_name>/merged_cf/combined_data/"$1"_combined_R1.fastq" \
"\nzcat "$1"*R2*fastq > ~/DATA/RNA-seq/<Project_name>/merged_cf/combined_data/"$1"_combined_R2.fastq" \
"\ncd ~/DATA/RNA-seq/<Project_name>/merged_cf/combined_data" \
"\n gzip "$1"*.fastq"}' >> zcat_reps.slurm
