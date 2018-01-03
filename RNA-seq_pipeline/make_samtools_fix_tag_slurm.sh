echo "#!/bin/bash
#SBATCH --job-name=samtools_fix_tag
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --partition=general
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=user@uchc.edu
#SBATCH -o samtools_fix_tag_%j.out
#SBATCH -e samtools_fix_tag_%j.err

source ~/.bashrc 
cd ~/ANALYSIS/RNA-seq/<Project_name>/merged_cf
export ANALYSISDIR=~/ANALYSIS/RNA-seq/<Project_name>/merged_cf"  > samtools_fix_tag.slurm

cat sample_list.txt | awk '{
print "\n# include reads that are 2nd in a pair (128);" \
"\n# exclude reads that are mapped to the reverse strand (16)" \
"\nsamtools view -b -f 128 -F 16 \$ANALYSISDIR/"$1"_PE/accepted_hits.bam > "$1".fwd1.bam" \
"\n" \
"\n# exclude reads that are mapped to the reverse strand (16) and" \
"\n# first in a pair (64): 64 + 16 = 80" \
"\nsamtools view -b -f 80 \$ANALYSISDIR/"$1"_PE/accepted_hits.bam > "$1".fwd2.bam" \
"\n" \
"\n# combine the temporary files" \
"\nsamtools merge -f "$1".fwd.bam "$1".fwd1.bam "$1".fwd2.bam" \
"\n" \
"\n# index the filtered BAM file" \
"\nsamtools index "$1".fwd.bam" \
"\n" \
"\n# run bamCoverage (only if DeepTools is able to run)" \
"\n#bamCoverage -b "$1".fwd.bam -o "$1".fwd.bigWig" \
"\n" \
"\n# remove the temporary files" \
"\n#rm "$1".fwd*.bam" \
"\n" \
"\n# if DeepTools cannot be run, use genomecov instead" \
"\nbedtools genomecov -ibam "$1".fwd.bam -bg -split > "$1"_pos.bg" \
"\nbedGraphToBigWig "$1"_pos.bg http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes 2017-10-24_"$1"_pos.bw" \
"\n" \
"\n#To get the file for transcripts that originated from the reverse strand:" \
"\n" \
"\n# include reads that map to the reverse strand (128)" \
"\n# and are second in a pair (16): 128 + 16 = 144" \
"\nsamtools view -b -f 144 \$ANALYSISDIR/"$1"_PE/accepted_hits.bam > "$1".rev1.bam" \
"\n" \
"\n# include reads that are first in a pair (64), but" \
"\n# exclude those ones that map to the reverse strand (16)" \
"\nsamtools view -b -f 64 -F 16 \$ANALYSISDIR/"$1"_PE/accepted_hits.bam > "$1".rev2.bam" \
"\n" \
"\n# merge the temporary files" \
"\nsamtools merge -f "$1".rev.bam "$1".rev1.bam "$1".rev2.bam" \
"\n" \
"\n# index the merged, filtered BAM file" \
"\nsamtools index "$1".rev.bam" \
"\n" \
"\n# run bamCoverage" \
"\n#bamCoverage -b "$1".rev.bam -o "$1".rev.bw" \
"\n" \
"\n# remove temporary files" \
"\n#rm "$1".rev*.bam" \
"\n" \
"\n# use genomecov instead" \
"\nbedtools genomecov -ibam "$1".rev.bam -bg -split > "$1"_neg.bg" \
"\nbedGraphToBigWig "$1"_neg.bg http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes "$1"_neg.bw"}' >> samtools_fix_tag.slurm
