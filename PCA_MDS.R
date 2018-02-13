library(SummarizedExperiment)
all_counts <- cbind(P6_counts,P13_counts,L1_counts)
all_counts <- as.data.frame(all_counts)
write.csv(all_counts, "/Users/belfordak/Desktop/Jan_seq_all/all_counts.csv")

setwd("/Users/belfordak/Desktop/Jan_seq_all")


colnames(all_counts) <- c("P6-r1_exon", "P6-r2_exon", "P6-r3_exon", "P6-r4_exon",
                          "P6-r1_intron", "P6-r2_intron", "P6-r3_intron", "P6-r4_intron",
                          "P13-r1_exon", "P13-r2_exon", "P13-r3_exon", "P13-r4_exon",
                          "P13-r1_intron", "P13-r2_intron", "P13-r3_intron", "P13-r4_intron",
                          "L1-r1_exon", "L1-r2_exon", "L1-r3_exon", "L1-r4_exon",
                          "L1-r1_intron", "L1-r2_intron", "L1-r3_intron", "L1-r4_intron")

allcounts_matrix <- as.matrix(all_counts)
coldata <- colData(all_counts)
countdata <- all_counts
coldata <- pData(all_counts)[,c("condition", "replicate")]

dds <- DESeqDataSet(allcounts_matrix, design = ~condition + replicate)

condition <- data.frame(row.names = c("P6-r1_exon", "P6-r2_exon", "P6-r3_exon", "P6-r4_exon",
                                    "P6-r1_intron", "P6-r2_intron", "P6-r3_intron", "P6-r4_intron",
                                    "P13-r1_exon", "P13-r2_exon", "P13-r3_exon", "P13-r4_exon",
                                    "P13-r1_intron", "P13-r2_intron", "P13-r3_intron", "P13-r4_intron",
                                    "L1-r1_exon", "L1-r2_exon", "L1-r3_exon", "L1-r4_exon",
                                    "L1-r1_intron", "L1-r2_intron", "L1-r3_intron", "L1-r4_intron"), condition=as.factor(c(rep("P6",4),rep("P13",4),rep("L1",4))))
                      
library(Biobase)
library(DEFormats)
library(DESeq2)
library(magrittr)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(ggrepel)

counts_as_RSE <- SummarizedExperiment(allcounts_matrix)                    
colData(counts_as_RSE)
countdata <- assay(counts_as_RSE)
head(countdata,3)
coldata <- colData(counts_as_RSE)
head(coldata)

condition=as.factor(c(rep("P6",8),rep("P13",8),rep("L1",8)))
mycols = data.frame(row.names = c("P6-r1_exon", "P6-r2_exon", "P6-r3_exon", "P6-r4_exon",
                                   "P6-r1_intron", "P6-r2_intron", "P6-r3_intron", "P6-r4_intron",
                                   "P13-r1_exon", "P13-r2_exon", "P13-r3_exon", "P13-r4_exon",
                                   "P13-r1_intron", "P13-r2_intron", "P13-r3_intron", "P13-r4_intron",
                                   "L1-r1_exon", "L1-r2_exon", "L1-r3_exon", "L1-r4_exon",
                                   "L1-r1_intron", "L1-r2_intron", "L1-r3_intron", "L1-r4_intron"), condition)

dds <- DESeqDataSetFromMatrix(countData = countdata, colData = mycols, design = ~condition)

rld <- rlog(dds, blind = FALSE)

colData(dds)

#for adding sample labels to plots
#+ geom_label_repel(aes(label =colnames(rld)), size = 1.8, label.size = 0, fill = "white")

#for PCA
plotPCA(rld)
rdl.sub.exon <- rld[,c(1:4,9:12,17:20)]
rdl.sub.intron <- rld[,c(5:8,13:16,21:24)]
head(rdl.sub.exon)
plotPCA(rdl.sub.exon) + geom_label_repel(aes(label =colnames(rdl.sub.exon)), size = 2.5, label.size = 0, fill = "white") + 
  ggtitle("Exon PCA with labels")
plotPCA(rdl.sub.intron) + geom_label_repel(aes(label =colnames(rdl.sub.intron)), size = 2.5, label.size = 0, fill = "white") + 
  ggtitle("Intron PCA with labels")

#for heatmap then mds
sampleDists <- dist(t(assay(rld)))

sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette(rev(brewer.pal(9, "Reds")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

mds <- as.data.frame(colData(rld))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color=condition)) +
  geom_point(size = 3, shape=3) + coord_fixed() + 
   ggtitle("Full MDS", subtitle = "exons on the left, introns on the right")

#exon mds
sampleDists.exon <- dist(t(assay(rdl.sub.exon)))
sampleDistMatrix.exon <- as.matrix(sampleDists.exon)
mds.exon <- as.data.frame(colData(rdl.sub.exon))  %>%
  cbind(cmdscale(sampleDistMatrix.exon))
ggplot(mds.exon, aes(x = `1`, y = `2`, color=condition)) +
  geom_point(size = 3, shape=3) + coord_fixed() + geom_label_repel(aes(label =colnames(rdl.sub.exon)), size = 2.5, label.size = 0, fill = "white") +
  ggtitle("exon MDS with labels")

#intron mds
sampleDists.intron <- dist(t(assay(rdl.sub.intron)))
sampleDistMatrix.intron <- as.matrix(sampleDists.intron)
mds.intron <- as.data.frame(colData(rdl.sub.intron))  %>%
  cbind(cmdscale(sampleDistMatrix.intron))
ggplot(mds.intron, aes(x = `1`, y = `2`, color=condition)) +
  geom_point(size = 3, shape=3) + coord_fixed() + geom_label_repel(aes(label =colnames(rdl.sub.intron)), size = 2.5, label.size = 0, fill = "white") + 
  ggtitle("intron MDS with labels")

