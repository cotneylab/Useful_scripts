# now in R
library(sva)
#library(pamr)
#library(limma)
#library(Biobase)
#library(BiocGenerics)
library(DESeq2)

# # Use DESeq2 to create a matrix, normalize and transform the counts


countsfile <- ("/home/CAM/awilderman/ANALYSIS/RNA-seq/2017-10-24_human_cf/combined_data/tophat_pipeline/R_data/reordered_CF_RNAseq_hg19_GCv10.txt")
countdata <- read.table(countsfile, header=TRUE, row.names=1)

# convert to matrix
countdata <- as.matrix(countdata)
head(countdata)

# assign condition (if putting these in this order, the columns in the counts data have to be in this order, too.)
condition <- factor(c(rep("CS13", 3), rep("CS15", 3), rep("CS14", 3), rep("CS17", 3)))

(coldata <- data.frame(row.names=colnames(countdata), condition))

# create the object
ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                  colData = coldata,
                                  design = ~ condition)
                                  
# transform for clustering (try three different methods and graph them to compare to eachother)
rld <- rlog(ddsMat, blind = FALSE)
head(assay(rld), 3)

vsd <- vst(ddsMat, blind = FALSE)
head(assay(vsd), 3)

# To plot how each transformation compares to the traditional Log2+1:
library("dplyr")
library("ggplot2")

ddsMat <- estimateSizeFactors(ddsMat)
# For the log2 approach, we need to first estimate size factors to account for sequencing depth, and then specify normalized=TRUE. Sequencing depth correction is done automatically for the rlog and the vst.

df <- bind_rows(
  as_data_frame(log2(counts(ddsMat, normalized=TRUE)[, 1:2]+1)) %>%
         mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))
  
colnames(df)[1:2] <- c("x", "y")  

pdf("CF_compare_transformation.pdf")                                
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation) 
dev.off()

# rld happens to come out the best, continue with this 

# Euclidean distance between samples to gauge similarity
sampleDists <- dist(t(assay(rld)))
sampleDists

# Make it a heatmap
library("pheatmap")
library("RColorBrewer")

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$condition )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf("CF_compare_sample_distances.pdf")                                
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
dev.off()

# Can also heatmap the Poisson Distance of the original (non-transformed)counts
library("PoiClaClu")
poisd <- PoissonDistance(t(counts(ddsMat)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( rld$condition )
colnames(samplePoisDistMatrix) <- NULL
pdf("CF_Poisson_distance_heatmap.pdf")
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)
dev.off()


# PCA plot
pdf("CF_PCA.pdf")                                
plotPCA(rld, intgroup = c("condition"))
dev.off()

# # Estimate surrogate variables

# the normalized data for sva (but what is this normalization?)
# and take normalized counts for which the average across samples is >1
# the known variable is condition (carnegie stage)
dat <- counts(ddsMat, normalized = TRUE)
idx <- rowMeans(dat) > 1
dat <- dat[idx,]
mod <- model.matrix (~ condition, colData(ddsMat))
mod0 <- model.matrix (~ 1, colData(ddsMat))
n.sv = num.sv(dat,mod,method="leek")
svseq <- svaseq(dat, mod, mod0)

# because 3 significant surrogate variables were estimated, plot these
pdf("CF_sva.pdf")
par(mfrow = c(2, 1), mar = c(3,5,3,1))
for (i in 1:3) {
  stripchart(svseq$sv[, i] ~ ddsMat$condition, vertical = TRUE, main = paste0("SV", i))
  abline(h = 0)
 }
dev.off()


ddssva <- ddsMat
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
ddssva$SV3 <- svseq$sv[,3]
design(ddssva) <- ~ SV1 + SV2 + SV3 + condition

ddssva <- DESeq(ddssva)
ressva <- results(ddssva)
ddssva
ressva 

ddsMat <- DESeq(ddsMat)
res <- results(ddsMat)
ddsMat
res 

# now we have two results sets, one corrected for surrogate variables (ressva) and one not (res), the default is to compare the first and last group, so CS17 vs CS13, but we can specify and name this comparison and others using contrast and give each set a unique and easy name

# Use both sets of results to do pairwise comparisons of gene expression 

# Annotating results will make it easier to search by gene name, load necessary packages first:
library("AnnotationDbi")
library("org.Hs.eg.db")

# to compare groups of interest use contrast:
# Comparisons of the data without surrogate variable correction (use ddsMat)
CS17_v_CS13_res <- results(ddsMat, contrast=c("condition","CS17","CS13")) 
table(CS17_v_CS13_res$padj<0.01)
# Annotate
CS17_v_CS13_res$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(CS17_v_CS13_res),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
CS17_v_CS13_res$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(CS17_v_CS13_res),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
## Order by adjusted p-value
CS17_v_CS13_res <- CS17_v_CS13_res[order(CS17_v_CS13_res$padj), ]
## Merge with normalized count data
CS17_v_CS13_resdata <- merge(as.data.frame(CS17_v_CS13_res), as.data.frame(counts(ddsMat, normalized=TRUE)), by="row.names", sort=FALSE)
names(CS17_v_CS13_resdata)[1] <- "Gene"
head(CS17_v_CS13_resdata)
## Write results
write.csv(CS17_v_CS13_resdata, file="CF_Tissue-CS17_v_CS13_diffexpr-results.csv") # table with log-2-fold change in expression of CF tissue stage over other, B-H adjusted p-value, ordered with lowest adjusted p-value on top (a list of everything and the DESeq2 results for that comparison)


# # Do this for other groups of interest, too:
CS15_v_CS13_res <- results(ddsMat, contrast=c("condition","CS15","CS13")) 
table(CS15_v_CS13_res$padj<0.01)
# Annotate
CS15_v_CS13_res$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(CS15_v_CS13_res),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
CS15_v_CS13_res$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(CS15_v_CS13_res),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
## Order by adjusted p-value
CS15_v_CS13_res <- CS15_v_CS13_res[order(CS15_v_CS13_res$padj), ]
## Merge with normalized count data
CS15_v_CS13_resdata <- merge(as.data.frame(CS15_v_CS13_res), as.data.frame(counts(ddsMat, normalized=TRUE)), by="row.names", sort=FALSE)
names(CS15_v_CS13_resdata)[1] <- "Gene"
head(CS15_v_CS13_resdata)
## Write results
write.csv(CS15_v_CS13_resdata, file="CF_Tissue-CS15_v_CS13_diffexpr-results.csv") 

# #
CS14_v_CS13_res <- results(ddsMat, contrast=c("condition","CS14","CS13")) 
table(CS14_v_CS13_res$padj<0.01)
# Annotate
CS14_v_CS13_res$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(CS14_v_CS13_res),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
CS14_v_CS13_res$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(CS14_v_CS13_res),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
## Order by adjusted p-value
CS14_v_CS13_res <- CS14_v_CS13_res[order(CS14_v_CS13_res$padj), ]
## Merge with normalized count data
CS14_v_CS13_resdata <- merge(as.data.frame(CS14_v_CS13_res), as.data.frame(counts(ddsMat, normalized=TRUE)), by="row.names", sort=FALSE)
names(CS14_v_CS13_resdata)[1] <- "Gene"
head(CS14_v_CS13_resdata)
## Write results
write.csv(CS14_v_CS13_resdata, file="CF_Tissue-CS14_v_CS13_diffexpr-results.csv") 

# # 
CS17_v_CS14_res <- results(ddsMat, contrast=c("condition","CS17","CS14")) 
table(CS17_v_CS14_res$padj<0.01)
# Annotate
CS17_v_CS14_res$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(CS17_v_CS14_res),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
CS17_v_CS14_res$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(CS17_v_CS14_res),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
## Order by adjusted p-value
CS17_v_CS14_res <- CS17_v_CS14_res[order(CS17_v_CS14_res$padj), ]
## Merge with normalized count data
CS17_v_CS14_resdata <- merge(as.data.frame(CS17_v_CS14_res), as.data.frame(counts(ddsMat, normalized=TRUE)), by="row.names", sort=FALSE)
names(CS17_v_CS14_resdata)[1] <- "Gene"
head(CS17_v_CS14_resdata)
## Write results
write.csv(CS17_v_CS14_resdata, file="CF_Tissue-CS17_v_CS14_diffexpr-results.csv") 

# #
CS17_v_CS15_res <- results(ddsMat, contrast=c("condition","CS17","CS15")) 
table(CS17_v_CS15_res$padj<0.01)
# Annotate
CS17_v_CS15_res$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(CS17_v_CS15_res),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
CS17_v_CS15_res$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(CS17_v_CS15_res),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
## Order by adjusted p-value
CS17_v_CS15_res <- CS17_v_CS15_res[order(CS17_v_CS15_res$padj), ]
## Merge with normalized count data
CS17_v_CS15_resdata <- merge(as.data.frame(CS17_v_CS15_res), as.data.frame(counts(ddsMat, normalized=TRUE)), by="row.names", sort=FALSE)
names(CS17_v_CS15_resdata)[1] <- "Gene"
head(CS17_v_CS15_resdata)
## Write results
write.csv(CS17_v_CS15_resdata, file="CF_Tissue-CS17_v_CS15_diffexpr-results.csv") 
# # 
CS15_v_CS14_res <- results(ddsMat, contrast=c("condition","CS15","CS14")) 
table(CS15_v_CS14_res$padj<0.01)
# Annotate
CS15_v_CS14_res$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(CS15_v_CS14_res),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
CS15_v_CS14_res$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(CS15_v_CS14_res),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
## Order by adjusted p-value
CS15_v_CS14_res <- CS15_v_CS14_res[order(CS15_v_CS14_res$padj), ]
## Merge with normalized count data
CS15_v_CS14_resdata <- merge(as.data.frame(CS15_v_CS14_res), as.data.frame(counts(ddsMat, normalized=TRUE)), by="row.names", sort=FALSE)
names(CS15_v_CS14_resdata)[1] <- "Gene"
head(CS15_v_CS14_resdata)
## Write results
write.csv(CS15_v_CS14_resdata, file="CF_Tissue-CS15_v_CS14_diffexpr-results.csv") 







# #
# Comparisons of the data without surrogate variable correction (use ddssva)
CS17_v_CS13_ressva <- results(ddssva, contrast=c("condition","CS17","CS13")) 
table(CS17_v_CS13_ressva$padj<0.01)
# Annotate
CS17_v_CS13_ressva$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(CS17_v_CS13_ressva),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
CS17_v_CS13_ressva$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(CS17_v_CS13_ressva),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
## Order by adjusted p-value
CS17_v_CS13_ressva <- CS17_v_CS13_ressva[order(CS17_v_CS13_ressva$padj), ]
## Merge with normalized count data
CS17_v_CS13_ressvadata <- merge(as.data.frame(CS17_v_CS13_ressva), as.data.frame(counts(ddssva, normalized=TRUE)), by="row.names", sort=FALSE)
names(CS17_v_CS13_ressvadata)[1] <- "Gene"
head(CS17_v_CS13_ressvadata)
## Write results
write.csv(CS17_v_CS13_ressvadata, file="sva_CF_Tissue-CS17_v_CS13_diffexpr-results.csv") 


# # Do this for other groups of interest, too:
CS15_v_CS13_ressva <- results(ddssva, contrast=c("condition","CS15","CS13")) 
table(CS15_v_CS13_ressva$padj<0.01)
# Annotate
CS15_v_CS13_ressva$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(CS15_v_CS13_ressva),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
CS15_v_CS13_ressva$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(CS15_v_CS13_ressva),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
## Order by adjusted p-value
CS15_v_CS13_ressva <- CS15_v_CS13_ressva[order(CS15_v_CS13_ressva$padj), ]
## Merge with normalized count data
CS15_v_CS13_ressvadata <- merge(as.data.frame(CS15_v_CS13_ressva), as.data.frame(counts(ddssva, normalized=TRUE)), by="row.names", sort=FALSE)
names(CS15_v_CS13_ressvadata)[1] <- "Gene"
head(CS15_v_CS13_ressvadata)
## Write results
write.csv(CS15_v_CS13_ressvadata, file="sva_CF_Tissue-CS15_v_CS13_diffexpr-results.csv") 
# #
CS14_v_CS13_ressva <- results(ddssva, contrast=c("condition","CS14","CS13")) 
table(CS14_v_CS13_ressva$padj<0.01)
# Annotate
CS14_v_CS13_ressva$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(CS14_v_CS13_ressva),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
CS14_v_CS13_ressva$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(CS14_v_CS13_ressva),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
## Order by adjusted p-value
CS14_v_CS13_ressva <- CS14_v_CS13_ressva[order(CS14_v_CS13_ressva$padj), ]
## Merge with normalized count data
CS14_v_CS13_ressvadata <- merge(as.data.frame(CS14_v_CS13_ressva), as.data.frame(counts(ddssva, normalized=TRUE)), by="row.names", sort=FALSE)
names(CS14_v_CS13_ressvadata)[1] <- "Gene"
head(CS14_v_CS13_ressvadata)
## Write results
write.csv(CS14_v_CS13_ressvadata, file="sva_CF_Tissue-CS14_v_CS13_diffexpr-results.csv") 

# # 
CS17_v_CS14_ressva <- results(ddssva, contrast=c("condition","CS17","CS14")) 
table(CS17_v_CS14_ressva$padj<0.01)
# Annotate
CS17_v_CS14_ressva$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(CS17_v_CS14_ressva),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
CS17_v_CS14_ressva$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(CS17_v_CS14_ressva),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
## Order by adjusted p-value
CS17_v_CS14_ressva <- CS17_v_CS14_ressva[order(CS17_v_CS14_ressva$padj), ]
## Merge with normalized count data
CS17_v_CS14_ressvadata <- merge(as.data.frame(CS17_v_CS14_ressva), as.data.frame(counts(ddssva, normalized=TRUE)), by="row.names", sort=FALSE)
names(CS17_v_CS14_ressvadata)[1] <- "Gene"
head(CS17_v_CS14_ressvadata)
## Write results
write.csv(CS17_v_CS14_ressvadata, file="sva_CF_Tissue-CS17_v_CS14_diffexpr-results.csv") 

# #
CS17_v_CS15_ressva <- results(ddssva, contrast=c("condition","CS17","CS15")) 
table(CS17_v_CS15_ressva$padj<0.01)
# Annotate
CS17_v_CS15_ressva$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(CS17_v_CS15_ressva),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
CS17_v_CS15_ressva$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(CS17_v_CS15_ressva),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
## Order by adjusted p-value
CS17_v_CS15_ressva <- CS17_v_CS15_ressva[order(CS17_v_CS15_ressva$padj), ]
## Merge with normalized count data
CS17_v_CS15_ressvadata <- merge(as.data.frame(CS17_v_CS15_ressva), as.data.frame(counts(ddssva, normalized=TRUE)), by="row.names", sort=FALSE)
names(CS17_v_CS15_ressvadata)[1] <- "Gene"
head(CS17_v_CS15_ressvadata)
## Write results
write.csv(CS17_v_CS15_ressvadata, file="sva_CF_Tissue-CS17_v_CS15_diffexpr-results.csv") 

# # 
CS15_v_CS14_ressva <- results(ddssva, contrast=c("condition","CS15","CS14")) 
table(CS15_v_CS14_ressva$padj<0.01)
# Annotate
CS15_v_CS14_ressva$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(CS15_v_CS14_ressva),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
CS15_v_CS14_ressva$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(CS15_v_CS14_ressva),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
## Order by adjusted p-value
CS15_v_CS14_ressva <- CS15_v_CS14_ressva[order(CS15_v_CS14_ressva$padj), ]
## Merge with normalized count data
CS15_v_CS14_ressvadata <- merge(as.data.frame(CS15_v_CS14_ressva), as.data.frame(counts(ddssva, normalized=TRUE)), by="row.names", sort=FALSE)
names(CS15_v_CS14_ressvadata)[1] <- "Gene"
head(CS15_v_CS14_ressvadata)
## Write results
write.csv(CS15_v_CS14_ressvadata, file="sva_CF_Tissue-CS15_v_CS14_diffexpr-results.csv") 


