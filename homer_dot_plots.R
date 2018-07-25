library(homerkit)
library(ggplot2)
library(ggrepel)
library(dplyr)
#load homer outputs
N2_UP_H3K27ac_enhancer <- read_homer_output("N2_UP_H3K27ac_enhancer_homer")
N4_UP_H3K27ac_enhancer <- read_homer_output("N4_UP_H3K27ac_enhancer_homer")
#extract 
N2_UP_H3K27ac_enhancer_known = as.data.frame(N2_UP_H3K27ac_enhancer$known_motif_table[,])
N4_UP_H3K27ac_enhancer_known = as.data.frame(N4_UP_H3K27ac_enhancer$known_motif_table[,])
combined <- merge(N2_UP_H3K27ac_enhancer_known, N4_UP_H3K27ac_enhancer_known, by="motif_name")
combined$motif_name<-substr(combined$motif_name,1,8)
p <- ggplot(combined, aes(log10(combined$number_of_target_sequences_with_motif_of_4902),(-log10(combined$p_value.x))-(-log10(combined$p_value.y))))
p2 <- p + geom_point(aes(col=(-log10(combined$p_value.x))-(-log10(combined$p_value.y)),size=number_of_target_sequences_with_motif_of_4902/4902)) +scale_color_gradient2(midpoint = 0, low="blue",mid="grey", high="red")+theme(axis.text.y=element_text (size=6)) + geom_text_repel(aes(log10(combined$number_of_target_sequences_with_motif_of_4902),(-log10(combined$p_value.x))-(-log10(combined$p_value.y)), label= ifelse((-log10(combined$p_value.x))-(-log10(combined$p_value.y)) > 3 | (-log10(combined$p_value.x))-(-log10(combined$p_value.y)) < -3 ,as.character(combined$motif_name),'')))
pdf(file="2N_vs_4N_H3K27ac_enhancer_motif_2.pdf", 11.5,8)
p2
dev.off()
#load homer outputs
N2_UP_H3K4me2_enhancer <- read_homer_output("N2_UP_H3K4me2_enhancer_homer")
N4_UP_H3K4me2_enhancer <- read_homer_output("N4_UP_H3K4me2_enhancer_homer")
#extract 
N2_UP_H3K4me2_enhancer_known = as.data.frame(N2_UP_H3K4me2_enhancer$known_motif_table[,])
N4_UP_H3K4me2_enhancer_known = as.data.frame(N4_UP_H3K4me2_enhancer$known_motif_table[,])
combined <- merge(N2_UP_H3K4me2_enhancer_known, N4_UP_H3K4me2_enhancer_known, by="motif_name")
combined$motif_name<-substr(combined$motif_name,1,8)
p <- ggplot(combined, aes(log10(combined$number_of_target_sequences_with_motif_of_3858),(-log10(combined$p_value.x))-(-log10(combined$p_value.y))))
p2 <- p + geom_point(aes(col=(-log10(combined$p_value.x))-(-log10(combined$p_value.y)),size=number_of_target_sequences_with_motif_of_3858/3858)) +scale_color_gradient2(midpoint = 0, low="blue",mid="grey", high="red")+theme(axis.text.y=element_text (size=6)) + geom_text_repel(aes(log10(combined$number_of_target_sequences_with_motif_of_3858),(-log10(combined$p_value.x))-(-log10(combined$p_value.y)), label= ifelse((-log10(combined$p_value.x))-(-log10(combined$p_value.y)) > 3 | (-log10(combined$p_value.x))-(-log10(combined$p_value.y)) < -3 ,as.character(combined$motif_name),'')))
pdf(file="2N_vs_4N_H3K4me2_enhancer_motif.pdf", 11.5,8)
p2
dev.off()
#load homer outputs
N2_UP_H3K4me3_enhancer <- read_homer_output("N2_UP_H3K4me3_enhancer_homer")
N4_UP_H3K4me3_enhancer <- read_homer_output("N4_UP_H3K4me3_enhancer_homer")
#extract 
N2_UP_H3K4me3_enhancer_known = as.data.frame(N2_UP_H3K4me3_enhancer$known_motif_table[,])
N4_UP_H3K4me3_enhancer_known = as.data.frame(N4_UP_H3K4me3_enhancer$known_motif_table[,])
combined <- merge(N2_UP_H3K4me3_enhancer_known, N4_UP_H3K4me3_enhancer_known, by="motif_name")
combined$motif_name<-substr(combined$motif_name,1,8)
p <- ggplot(combined, aes(log10(combined$number_of_target_sequences_with_motif_of_98),(-log10(combined$p_value.x))-(-log10(combined$p_value.y))))
p2 <- p + geom_point(aes(col=(-log10(combined$p_value.x))-(-log10(combined$p_value.y)),size=number_of_target_sequences_with_motif_of_98/98)) +scale_color_gradient2(midpoint = 0, low="blue",mid="grey", high="red")+theme(axis.text.y=element_text (size=6)) + geom_text_repel(aes(log10(combined$number_of_target_sequences_with_motif_of_98),(-log10(combined$p_value.x))-(-log10(combined$p_value.y)), label= ifelse((-log10(combined$p_value.x))-(-log10(combined$p_value.y)) >= 2 | (-log10(combined$p_value.x))-(-log10(combined$p_value.y)) <= -2 ,as.character(combined$motif_name),'')))
pdf(file="2N_vs_4N_H3K4me3_enhancer_motif.pdf", 11.5,8)
p2
dev.off()
#load homer outputs
N2_UP_H3K27ac_promoter <- read_homer_output("N2_UP_H3K27ac_promoter_homer")
N4_UP_H3K27ac_promoter <- read_homer_output("N4_UP_H3K27ac_promoter_homer")
#extract 
N2_UP_H3K27ac_promoter_known = as.data.frame(N2_UP_H3K27ac_promoter$known_motif_table[,])
N4_UP_H3K27ac_promoter_known = as.data.frame(N4_UP_H3K27ac_promoter$known_motif_table[,])
combined <- merge(N2_UP_H3K27ac_promoter_known, N4_UP_H3K27ac_promoter_known, by="motif_name")
combined$motif_name<-substr(combined$motif_name,1,8)
p <- ggplot(combined, aes(log10(combined$number_of_target_sequences_with_motif_of_197),(-log10(combined$p_value.x))-(-log10(combined$p_value.y))))
p2 <- p + geom_point(aes(col=(-log10(combined$p_value.x))-(-log10(combined$p_value.y)),size=number_of_target_sequences_with_motif_of_197/197)) +scale_color_gradient2(midpoint = 0, low="blue",mid="grey", high="red")+theme(axis.text.y=element_text (size=6)) + geom_text_repel(aes(log10(combined$number_of_target_sequences_with_motif_of_197),(-log10(combined$p_value.x))-(-log10(combined$p_value.y)), label= ifelse((-log10(combined$p_value.x))-(-log10(combined$p_value.y)) >= 2 | (-log10(combined$p_value.x))-(-log10(combined$p_value.y)) <= -2 ,as.character(combined$motif_name),'')))
pdf(file="2N_vs_4N_H3K27ac_promoter_motif.pdf", 11.5,8)
p2
dev.off()
#load homer outputs
N2_UP_H3K4me2_promoter <- read_homer_output("N2_UP_H3K4me2_promoter_homer")
N4_UP_H3K4me2_promoter <- read_homer_output("N4_UP_H3K4me2_promoter_homer")
#extract 
N2_UP_H3K4me2_promoter_known = as.data.frame(N2_UP_H3K4me2_promoter$known_motif_table[,])
N4_UP_H3K4me2_promoter_known = as.data.frame(N4_UP_H3K4me2_promoter$known_motif_table[,])
combined <- merge(N2_UP_H3K4me2_promoter_known, N4_UP_H3K4me2_promoter_known, by="motif_name")
combined$motif_name<-substr(combined$motif_name,1,8)
p <- ggplot(combined, aes(log10(combined$number_of_target_sequences_with_motif_of_1116),(-log10(combined$p_value.x))-(-log10(combined$p_value.y))))
p2 <- p + geom_point(aes(col=(-log10(combined$p_value.x))-(-log10(combined$p_value.y)),size=number_of_target_sequences_with_motif_of_1116/1116)) +scale_color_gradient2(midpoint = 0, low="blue",mid="grey", high="red")+theme(axis.text.y=element_text (size=6)) + geom_text_repel(aes(log10(combined$number_of_target_sequences_with_motif_of_1116),(-log10(combined$p_value.x))-(-log10(combined$p_value.y)), label= ifelse((-log10(combined$p_value.x))-(-log10(combined$p_value.y)) >= 2 | (-log10(combined$p_value.x))-(-log10(combined$p_value.y)) <= -2 ,as.character(combined$motif_name),'')))
pdf(file="2N_vs_4N_H3K4me2_promoter_motif.pdf", 11.5,8)
p2
dev.off()
#load homer outputs
N2_UP_H3K4me3_promoter <- read_homer_output("N2_UP_H3K4me3_promoter_homer")
N4_UP_H3K4me3_promoter <- read_homer_output("N4_UP_H3K4me3_promoter_homer")
#extract 
N2_UP_H3K4me3_promoter_known = as.data.frame(N2_UP_H3K4me3_promoter$known_motif_table[,])
N4_UP_H3K4me3_promoter_known = as.data.frame(N4_UP_H3K4me3_promoter$known_motif_table[,])
combined <- merge(N2_UP_H3K4me3_promoter_known, N4_UP_H3K4me3_promoter_known, by="motif_name")
combined$motif_name<-substr(combined$motif_name,1,8)
p <- ggplot(combined, aes(log10(combined$number_of_target_sequences_with_motif_of_970),(-log10(combined$p_value.x))-(-log10(combined$p_value.y))))
p2 <- p + geom_point(aes(col=(-log10(combined$p_value.x))-(-log10(combined$p_value.y)),size=number_of_target_sequences_with_motif_of_970/970)) +scale_color_gradient2(midpoint = 0, low="blue",mid="grey", high="red")+theme(axis.text.y=element_text (size=6)) + geom_text_repel(aes(log10(combined$number_of_target_sequences_with_motif_of_970),(-log10(combined$p_value.x))-(-log10(combined$p_value.y)), label= ifelse((-log10(combined$p_value.x))-(-log10(combined$p_value.y)) >= 2 | (-log10(combined$p_value.x))-(-log10(combined$p_value.y)) <= -2 ,as.character(combined$motif_name),'')))
pdf(file="2N_vs_4N_H3K4me3_promoter_motif.pdf", 11.5,8)
p2
dev.off()