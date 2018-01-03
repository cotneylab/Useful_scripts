### make_featureCounts_hg19_GCv10_2_fr_slurm.sh
# can't run the make script until the sams exist!!
# this version uses hg19, Gencode GCv10 and assumes a library from the illumina Truseq stranded kit (2 fr)

export ANALYSISDIR=~/ANALYSIS/RNA-seq/<Project_name>/merged_cf
export GTFDIR=~/TOOLS/bcbio/genomes/Hsapiens/hg19/annotation

echo -e "#!/bin/bash
#SBATCH --job-name=feature_counts
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH --partition=general
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=user@uchc.edu
#SBATCH -o feature_counts_%j.out
#SBATCH -e feature_counts_%j.err
source ~/.bashrc
cd $ANALYSISDIR
featureCounts -s 2 -p -B -g gene_id -t exon -a $GTFDIR/gencode.v10.annotation.gtf -o counts_RNAseq_hg19_GCv10_2_fr.txt -C \\"}' > featureCounts_hg19_GCv10_2_fr.slurm
ls -1 *.sam | awk '{print $1 }' ORS=' ' >> featureCounts_hg19_GCv10_2_fr.slurm
# this last bit removes newline for list of sam files

# makes a single file counts_RNAseq_hg19_GCv10_2_fr.txt