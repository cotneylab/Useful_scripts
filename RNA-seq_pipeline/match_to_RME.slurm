#!/bin/bash
#SBATCH --job-name=match_to_RME
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH --partition=general
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=user@uchc.edu
#SBATCH -o match_to_RME_%j.out
#SBATCH -e match_to_RME_%j.err
source ~/.bashrc
cd ~/ANALYSIS/RNA-seq/<Project_name>/merged_cf

# Match Geneids between RNA-seq data and Roadmap Epigenomics files
# First sort both on same column
# To sort RNA-seq data and preserve header rows:
# use this one: retains header
(head -n 1 counts_RNAseq_hg19_GCv10_2_fr.txt && tail -n +3 counts_RNAseq_hg19_GCv10_2_fr.txt | sort) > sorted_counts_RNAseq_hg19_GCv10_2_fr.txt

# and sort Roadmap data, which only has one header line 
# This actually has to be the counts table, so download that and start from here:
# file name is 57epigenomes.N.pc
cp ~/DATA/RNA-seq/ROADMAP/57epigenomes.N.pc ~/ANALYSIS/<Project_name>/merged_cf/
cd ~/ANALYSIS/<Project_name>/merged_cf
(head -n 1 57epigenomes.N.pc && tail -n +2 57epigenomes.N.pc | sort) > sorted_57epigenomes.N.pc

# For join may want to specify a few things:
# first column is default, ok
# so many more are unpaired, this is probably a problem with the Geneid versions (. vs no .)
# This will change the Geneid to non-decimal, like in Roadmap and removes the intervening columns with chromosome number, start, end, strand and length, leaving only Geneid and reads along with the header column: 
cut -f1 sorted_counts_RNAseq_hg19_GCv10_2_fr.txt | sed 's/\./\t/g' | cut -f1 | paste - sorted_counts_RNAseq_hg19_GCv10_2_fr.txt | cut -f 1,8-19 > nodot_sorted_counts_RNAseq_hg19_GCv10_2_fr.txt
join -i --check-order --header -1 1 sorted_57epigenomes.N.pc nodot_sorted_counts_RNAseq_hg19_GCv10_2_fr.txt > new_joined_57epi_counts_RNAseq_hg19_GCv10_2rf.counts
join -v 2 -i --check-order --header -1 1 sorted_57epigenomes.N.pc nodot_sorted_counts_RNAseq_hg19_GCv10_2_fr.txt > new_unpaired_57epi_counts_RNAseq_hg19_GCv10_2rf.counts

