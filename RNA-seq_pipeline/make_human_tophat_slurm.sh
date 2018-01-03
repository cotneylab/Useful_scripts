echo "#!/bin/bash
#SBATCH --job-name=human_tophat
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --partition=general
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=user@uchc.edu
#SBATCH -o human_tophat_%j.out
#SBATCH -e human_tophat_%j.err

source ~/.bashrc 
cd ~/ANALYSIS/RNA-seq/<Project_name>/merged_cf
export ANALYSISDIR=~/ANALYSIS/RNA-seq/<Project_name>/merged_cf
export GTFDIR=~/TOOLS/bcbio/genomes/Hsapiens/hg19/annotation
export INDEXDIR=~/TOOLS/bcbio/genomes/Hsapiens/hg19/bowtie2
export FQDIR=~/ANALYSIS/RNA-seq/<Project_name>/merged_cf" > human_tophat.slurm

cat sample_list.txt | awk '{ \
print "\ntophat -o \$ANALYSISDIR/"$1"_PE --library-type fr-firststrand -p 16 -G \$GTFDIR/gencode.v10.annotation.gtf --no-discordant --mate-inner-dist 100 \$INDEXDIR/hg19 \$FQDIR/"$1"_forward_paired.fq.gz \$FQDIR/"$1"_reverse_paired.fq.gz" \
}' >> human_tophat.slurm
