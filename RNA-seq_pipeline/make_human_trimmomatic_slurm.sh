echo "#!/bin/bash
#SBATCH --job-name=human_trimmomatic
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --partition=general
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=user@uchc.edu
#SBATCH -o human_trimmomatic_%j.out
#SBATCH -e human_trimmomatic_%j.err
source ~/.bashrc
cd ~/ANALYSIS/RNA-seq/<Project_name>/merged_cf
export TRIMDIR=/isg/shared/apps/Trimmomatic/0.36
export ADAPTERDIR=/isg/shared/apps/Trimmomatic/0.36/adapters
export FASTQDIR=~/DATA/RNA-seq/<Project_name>/merged_cf
export ANALYSISDIR=~/ANALYSIS/RNA-seq/<Project_name>/merged_cf
module load Trimmomatic/0.36" > human_trimmomatic.slurm
cat sample_list.txt | awk '{ \
print "\njava -jar \$TRIMDIR/trimmomatic-0.36.jar PE -phred33 \$FASTQDIR/"$1"_R1.fastq.gz \$FASTQDIR/"$1"_R2.fastq.gz \$ANALYSISDIR/"$1"_forward_paired.fq.gz \$ANALYSISDIR/"$1"_forward_unpaired.fq.gz \$ANALYSISDIR/"$1"_reverse_paired.fq.gz $ANALYSISDIR/"$1"_reverse_unpaired.fq.gz ILLUMINACLIP:\$ADAPTERDIR/TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36" \
}' >> human_trimmomatic.slurm
