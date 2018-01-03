cat namelist.txt | awk '{print $1"_set60 <- read.table(\x22"$1"_bg_list.txt\x22, header=TRUE)" \
"\nhead("$1"_set60)" \
"\npdf(\x22"$3"_set60_Volcano_plot.pdf\x22)" \
"\nwith("$1"_set60, plot(log2FoldChange, -log10(padj), pch=20, main=\x22"$3" set60 Volcano plot\x22, xlim=c(-16,16)))" \
"\n" \
"\n# Add colored points: red if padj<0.01, orange of log2FC>1.5, green if both)" \
"\nwith(subset("$1"_set60, padj<.01 ), points(log2FoldChange, -log10(padj), pch=20, col=\x22red\x22))" \
"\nwith(subset("$1"_set60, abs(log2FoldChange)>1.5), points(log2FoldChange, -log10(padj), pch=20, col=\x22orange\x22))" \
"\nwith(subset("$1"_set60, padj<.01 & abs(log2FoldChange)>1.5), points(log2FoldChange, -log10(padj), pch=20, col=\x22green\x22))" \
"\n" \
"\n# Label points with the textxy function from the calibrate plot" \
"\nlibrary(calibrate)" \
"\nwith(subset("$1"_set60, padj<1.0e-5 & abs(log2FoldChange)>2), textxy(log2FoldChange, -log10(padj), labs=symbol, cex=.8))" \
"\ndev.off()" \
"\n# Out of the "set60" results, create a list of the upregulated genes" \
"\n"$2"_set60_Sig <- subset("$1"_set60, padj < 0.01)" \
"\nhead("$2"_set60_Sig[ order("$2"_set60_Sig$log2FoldChange, decreasing = TRUE), ])" \
"\n"$2"_set60_Sig_sorted <- (as.data.frame("$2"_set60_Sig[ order("$2"_set60_Sig$log2FoldChange, decreasing = TRUE), ]))" \
"\nwrite.csv("$2"_set60_Sig_sorted, file=\x22"$2"_set_60_Sig_sorted.csv\x22)" \
"\n" \
"\n"$2"_set60_Sig_upreg_in_"$4" <- subset("$2"_set60_Sig_sorted, log2FoldChange > 1)" \
"\nwrite.csv("$2"_set60_Sig_upreg_in_"$4", file=\x22"$2"_set60_Sig_upreg_in_"$4".csv\x22)" \
"\n"$2"_set60_Sig_downreg_in_"$4" <- subset("$2"_set60_Sig_sorted, log2FoldChange < -1)" \
"\nwrite.csv("$2"_set60_Sig_downreg_in_"$4", file=\x22"$2"_set60_Sig_downreg_in_"$4".csv\x22)" \
"\n" }' > CF_compare_volcano_with_set60_lists.R
