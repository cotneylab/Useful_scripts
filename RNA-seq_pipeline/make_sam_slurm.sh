export ANALYSISDIR=~/ANALYSIS/RNA-seq/<Project_name>/merged_cf

echo -e "#!/bin/bash
#SBATCH --job-name=make_sam
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH --partition=general
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=user@uchc.edu
#SBATCH -o make_sam_%j.out
#SBATCH -e make_sam_%j.err
source ~/.bashrc
cd $ANALYSISDIR" > make_sam.slurm

cat sample_list.txt | awk '{
print "\nsamtools sort -n -O SAM -o ~/ANALYSIS/RNA-seq/<Project_name>/merged_cf/"$1".sorted.sam ~/ANALYSIS/RNA-seq/<Project_name>/merged_cf/"$1"_accepted_hits.bam" }' >> make_sam.slurm

