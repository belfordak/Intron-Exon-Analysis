library(DESeq2)
library(GenomicAlignments)
library(GenomicFeatures)

#function to get the intervals of Grange object
interval <- function(gr){
  require(GenomicRanges)
  require(GenomicFeatures)
  len.gr <- length(gr)
  if (len.gr <= 1 ) {
    intron <- GRanges()}
  else {
    gap <- gaps(gr)
    intron <- gap[2:length(gap)]
  }
  return(intron)
}

#function to get the only exonic ranges, 
#query is a disjoined Grange object of a gene, exByTxUnlist is the GRange object constains all the exons of transcripts
exonOnlyRange <- function(query, exByTxUnlist){
  exonCounts <- countOverlaps(query, exByTxUnlist)
  return (query[exonCounts == max(exonCounts)])
}

#GrangeList to bedfile
grlToBed <- function(grl){
  gr <- unlist(grl)
  df <- data.frame(seqnames=seqnames(gr),
                   starts=start(gr)-1,
                   ends=end(gr),
                   names=names(gr),
                   scores=c(rep(".", length(gr))),
                   strands=strand(gr))
  
  write.table(df, file=paste0(substitute(grl),".bed"), quote=F, sep="\t", row.names=F, col.names=F)  
}

#GrangeList to GTF
grlToGTF <- function (grl,feature=c("exon","transcript","gene","intron")){
  gr <- unlist(grl)
  gene_id <- paste0('"',names(gr),'";')
  df <- data.frame(seqnames=seqnames(gr),
                   source=c(rep("NIDDK_LGP", length(gr))),
                   feature=feature,
                   starts=start(gr)-1,
                   ends=end(gr),
                   scores=c(rep(".", length(gr))),
                   strands=strand(gr),
                   frame=c(rep(".", length(gr))),
                   attribute=paste("gene_id",gene_id, "transcript_id", gene_id, "gene_name", gene_id, sep = " "))
  write.table(df, file=paste0(substitute(grl),".gtf"), quote=F, sep="\t", row.names = F, col.names = F)
}
#awk -F "\t" '$3 == "exon" { print $9}' Mus_musculus.GRCm38.90.gtf | tr -d ";\"" | awk -F " " '$16 == "protein_coding" {print $2}' | uniq > protein.coding.geneid
#selecting for only "protein coding" elimintes the calling of incorrect genes 
txdb <- makeTxDbFromGFF("/Users/belfordak/Data/Genome_files/UCSC_genes_mm10.gtf")
PC.geneid <- read.table("/Users/belfordak/Data/Genome_files/protein.coding.geneid", sep = '\t', header = F)$V1
library(mygene)
PC.genename <- queryMany(PC.geneid, scopes = "ensembl.gene", fields = "symbol")
PC.symbol <- unique(PC.genename$symbol)

#find out all the protein coding genes
exondb <- exonsBy(txdb, by= "gene", use.names=F) #GRangeList object of genes with overlapping exons
exondb.gr <- unlist(exondb)
exondb.pc <- subset(exondb.gr, names(exondb.gr) %in% PC.symbol)
exondb.pc <- split(exondb.pc, names(exondb.pc))
grlToBed(exondb.pc)

#filter out genes mapped to multiple places
exondb.filtered <- GRangesList()
for (i in 1:length(exondb.pc)) {
  if (length(unique(as.vector(strand(exondb.pc[[i]])))) > 1) next
  if (length(unique(as.vector(seqnames(exondb.pc[[i]])))) > 1) next
  exondb.filtered <- c(exondb.filtered, exondb.pc[i])
}


exonCollapsedb.filtered <- GenomicRanges::reduce(exondb.filtered) #GRangeList object of genes with non-overlapping exons
grlToGTF(exonCollapsedb.filtered,"exon")
exonCollapseTxdb <-makeTxDbFromGFF("exonCollapsedb.filtered.gtf")
introndb.filtered <- intronsByTranscript(exonCollapseTxdb, use.names=T)


names(introndb.filtered) <- paste0("Intron_",names(introndb.filtered))
exonPlusIntrondb.filtered <- c(exonCollapsedb.filtered, introndb.filtered) ##Here's the file for better library normalisation (accounts for both intron and exon, and within the context of the gene, for each file (intron and exon file))

##Here's the point that differs between different samples/conditions. The above code only need be run once (per enviornment)

setwd("/Users/belfordak/Desktop/Jan_seq_all/bams_bais/P6")

bmfls <- list.files(pattern="bam$",full.names=T)
bmfls_list <- BamFileList(bmfls)
geneCnt <- summarizeOverlaps(exonPlusIntrondb.filtered,bmfls_list,ignore.strand=T,singleEnd=TRUE,fragments=FALSE)
samples <- colnames(geneCnt)
condition <- c("P6_r1", "P6_r2", "P6_r3","P6_r4")
colData(geneCnt) <- DataFrame(samples, condition, row.names = samples)
dds <- DESeqDataSet(geneCnt, design = ~condition)
dds <- DESeq(dds)
fpkm.all <- fpkm(dds)
intron.names <- grep("Intron_", row.names(fpkm.all))
exon.names <- grep("Intron_", row.names(fpkm.all), invert = T)
fpkm.intron <- fpkm.all[intron.names,]
fpkm.exon <- fpkm.all[exon.names,]
rownames(fpkm.intron) <- gsub("Intron_","",rownames(fpkm.intron))
fpkm.transform <- merge(fpkm.exon, fpkm.intron, by=0)
names(fpkm.transform) <- c("Gene", "P6_r1_exon_fpkm","P6_r2_exon_fpkm", "P6_r3_exon_fpkm", "P6_r4_exon_fpkm", 
                           "P6_r1_intron_fpkm", "P6_r2_intron_fpkm", "P6_r3_intron_fpkm", "P6_r4_intron_fpkm")

write.csv(fpkm.transform, file = "P6_lowseq_exon_intron_fpkm.csv", quote = F, row.names = F)

#for counts
library(GenomicAlignments)
library(GenomicRanges)
library(GenomicFeatures)
library(Rsamtools)
library(pbapply)
library(DESeq2)
library(DEFormats)

#Creating count values
intronCnt <- summarizeOverlaps(introndb.filtered, bmfls_list,ignore.strand=T,singleEnd=TRUE,fragments=FALSE)
exonCnt <- summarizeOverlaps(exonCollapsedb.filtered, bmfls_list,ignore.strand=T,singleEnd=TRUE,fragments=FALSE)

#Creating intron count table
colData(intronCnt)$condition = condition
ok <- rep(c("1", "2", "3", "4"))
colData(intronCnt)$condition = ok
dds_intron_datset = DESeqDataSet(intronCnt, design = ~ condition)
dge_intron_lst = DGEList(dds_intron_datset)
dge_intron_lst = as.DGEList(dds_intron_datset)
head(dge_intron_lst)

P6_intron_table <- as.matrix(dge_intron_lst)

#Creating exon count table
colData(exonCnt)$condition = ok
dds_exon_datset = DESeqDataSet(exonCnt, design = ~ condition)
dge_exon_lst = DGEList(dds_exon_datset)
ok <- rep(c("1", "2", "3", "4"))
colData(exonCnt)$condition = ok
dds_exon_datset = DESeqDataSet(exonCnt, design = ~ condition)
dge_exon_lst = as.DGEList(dds_exon_datset)


P6_exon_table <- as.matrix(dge_exon_lst)
P6_counts <- cbind(P6_exon_table,P6_intron_table)
P6_counts_df <- as.data.frame(P6_counts)
head(P13_counts)
write.csv(P6_counts_df, file = "P6_counts.csv")

