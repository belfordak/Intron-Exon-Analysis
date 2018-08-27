library(magrittr)
library(dplyr)
library(DESeq2)
library(DEFormats)
library(RColorBrewer)
library(ggplot2)

setwd('/Users/belfordak/Desktop/L1:P6:P13_documents/Tcells_bams/')

base_csv_Tcell <- read.csv("new-lib-Tcell_exon_intron_fpkm.csv")
View(base_csv_Tcell)

colnames(base_csv_Tcell)[15:20] <- c("ratio_0hr-r1", "ratio_0hr-r2", "ratio_2w-r1", "ratio_2w-r2", "ratio_72hr-r1", "ratio_72hr-r2")

#Avg replicates

base_csv_Tcell$Mean_0hr_exon_fpkm <- rowMeans(base_csv_Tcell[,2:3])
base_csv_Tcell$Mean_2w_exon_fpkm <- rowMeans(base_csv_Tcell[,4:5])
base_csv_Tcell$Mean_72hr_exon_fpkm <- rowMeans(base_csv_Tcell[,6:7])
base_csv_Tcell$Mean_0hr_intron_fpkm <- rowMeans(base_csv_Tcell[,8:9])
base_csv_Tcell$Mean_2w_intron_fpkm <- rowMeans(base_csv_Tcell[,10:11])
base_csv_Tcell$Mean_72hr_intron_fpkm <- rowMeans(base_csv_Tcell[,12:13])

base_csv_Tcell <- mutate(base_csv_Tcell, avg_ratio_0hr = Mean_0hr_exon_fpkm / Mean_0hr_intron_fpkm, avg_ratio_2w = Mean_2w_exon_fpkm / Mean_2w_intron_fpkm, avg_ratio_72hr = Mean_72hr_exon_fpkm / Mean_72hr_intron_fpkm)

# Round 
base_csv_Tcell$avg_ratio_0hr<-round(base_csv_Tcell$avg_ratio_0hr, 3) 
base_csv_Tcell$avg_ratio_2w<-round(base_csv_Tcell$avg_ratio_2w, 3)
base_csv_Tcell$avg_ratio_72hr<-round(base_csv_Tcell$avg_ratio_72hr, 3)

View(base_csv_Tcell)

#Did this for excel -> R simplicity
base_csv_Tcell[base_csv_Tcell == '#VALUE!'] <- NA
base_csv_Tcell[base_csv_Tcell == '#DIV/0!'] <- NA
base_csv_Tcell[,14] <- NULL


#1. Apply min fpkm value of 2 (cuts out a lot of nonsignificant genes)
base_csv_Tcell <- as.data.frame(base_csv_Tcell)
trim_test <- base_csv_Tcell[rowSums(base_csv_Tcell[,2:7] < 2) <=1 , , drop = FALSE] #applies min cutoff for atleast one cell of exons to be 2, otherwise delete row (at least one replicate, of one sample must be >= 2)
trim_test2 <- trim_test[apply(trim_test, 1, function(y) !all(is.na(y))),] #would be for removing row if all are NAs. This doesn't work for us b/c of the above function takes care of it, but it's another option if you don't want the min fpkm filter applied
trim_3 <- trim_test[apply(trim_test[,14:19], 1, function(y) !all(is.na(y))),] #removes rows where all ratio counts are NAs

trim_4 <- trim_test[!is.na(trim_test[,14:19]$y),]

trim_4 <- trim_3[!is.na(trim_3[,14:19])]

trim_4 <- trim_3[complete.cases(trim_3),]                      
   
trim_4.5 <- trim_test[complete.cases(trim_test),] #just checking here

write.csv(trim_4, "~/Desktop/Tcell_only/noNA_fpkm.csv")

base_csv_Tcell <- read.table( "~/Desktop/Tcell_only/Tcell_master_p1.txt", header = TRUE)
base_csv_Tcell2 <- base_csv_Tcell[,-1]
rownames(base_csv_Tcell2) <- base_csv_Tcell[,1]
base_csv_Tcell <- base_csv_Tcell2


#Avg replicates

base_csv_Tcell$Mean_0hr_exon_fpkm <- rowMeans(base_csv_Tcell[,1:2])
base_csv_Tcell$Mean_2w_exon_fpkm <- rowMeans(base_csv_Tcell[,3:4])
base_csv_Tcell$Mean_72hr_exon_fpkm <- rowMeans(base_csv_Tcell[,5:6])
base_csv_Tcell$Mean_0hr_intron_fpkm <- rowMeans(base_csv_Tcell[,7:8])
base_csv_Tcell$Mean_2w_intron_fpkm <- rowMeans(base_csv_Tcell[,9:10])
base_csv_Tcell$Mean_72hr_intron_fpkm <- rowMeans(base_csv_Tcell[,11:12])

base_csv_Tcell$Mean_0hr_ratio <- rowMeans(base_csv_Tcell[,13:14])
base_csv_Tcell$Mean_2w_ratio <- rowMeans(base_csv_Tcell[,15:16])
base_csv_Tcell$Mean_72hr_ratio <- rowMeans(base_csv_Tcell[,17:18])

#Log transform

log_mean_basecsv_Tcell <- log2(base_csv_Tcell) 
log_mean_basecsv_Tcell <- round(log_mean_basecsv_Tcell, 3)
log_mean_basecsv_Tcell <- cbind(base_csv_Tcell[,1], log_mean_basecsv_Tcell)
colnames(base_csv_Tcell)[1] <- c("Gene")
View(log_mean_basecsv_Tcell)

#2. Subset data into ratio values of 1,5,10,100 - These are just notes to myself thinking about how to do this - I ended up doing it downstream though in excel
#newdat1.0 <- subset(base_csvvvv, column/variable? > 1.0, select=c(coluns to keep))
#newdat5.0 <- subset(base_csvvvv, column/variable? > 5.0, select=c(coluns to keep))
#newdat10.0 <- subset(base_csvvvv, column/variable? > 10.0, select=c(coluns to keep))
#newdat100.0 <- subset(base_csvvvv, column/variable? > 100.0, select=c(coluns to keep))

#GSEA for enrcihemnt of genes between conditions is good...random comment I guess?...but true...I did do this...why did I put this here?...but I guess it's going to stay...

Gene <- row.names(base_csv_Tcell)
Genes <- row.names(log_mean_basecsv_Tcell)

#4. Plots
#(I've kept in the "#old" code for personal reminders on older figures, but the "#NEW" and "#with log transform" are better and used in final figs)

#a) ratio vs. exon fpkm

#old
ggplot(basecsv_Tcell, aes(avg_ratio_0hr,Mean_0hr_exon_fpkm)) + geom_point(alpha=0.5) + 
  ggtitle("T cell exon/intron ratio to exon fpkm", subtitle = "applied min cutoff 2, not log transformed, replicates combined")

#NEW
ggplot(base_csv_Tcell, aes(Mean_72hr_ratio,Mean_72hr_exon_fpkm, label=Gene)) + geom_point(alpha=0.5) + 
  geom_text(aes(label=ifelse(Mean_72hr_exon_fpkm>700,as.character(Gene),'')),hjust=0,vjust=0) + 
  ggtitle("T cell 72 hr exon/intron ratio to exon fpkm", subtitle = "applied min exon cutoff 2, not log transformed, replicates combined")

#with log transform
ggplot(log_mean_basecsv_Tcell, aes(Mean_72hr_ratio, Mean_72hr_exon_fpkm, label=Genes)) + geom_point(alpha=0.5) + 
  geom_text_repel(aes(label=ifelse(Mean_72hr_exon_fpkm>9.5,as.character(Genes),'')),hjust=0,vjust=0) + 
  ggtitle("T cell 72 hr exon/intron ratio to exon fpkm", subtitle = "applied min exon cutoff 2, log transformed, replicates combined") +
  

#correlation factor
cor(log_mean_basecsv_Tcell$Mean_72hr_ratio, log_mean_basecsv_Tcell$Mean_72hr_exon_fpkm)

#old, useful, but not necessary with new, smaller data set
library(ggrepel)
ggplot(log_mean_base_csv_Tcell, aes(avg_ratio_0hr,Mean_0hr_exon_fpkm, label=Gene)) + geom_point(alpha=0.5) + 
  geom_text_repel(aes(label=ifelse(Mean_0hr_exon_fpkm>8,as.character(Gene),'')),hjust=0,vjust=0) + 
  ggtitle("T cell exon/intron ratio to exon fpkm", subtitle = "applied min cutoff X, log transformed, replicates combined")

#b) ratio vs exon fc
  #need to make x-axis = FC = ((B-A)/A) = (72-0)/0) = ((avg exon fpkm 72 - avg exon fpk 0) / (avg exon fpkm 0)...will need to create a variable..just new column for this..shouldn't be too hard... 

base_csv_Tcell$FC_72hr_exon <- ((base_csv_Tcell$Mean_72hr_exon_fpkm - base_csv_Tcell$Mean_0hr_exon_fpkm) / base_csv_Tcell$Mean_0hr_exon_fpkm)
#options?
#log_mean_basecsv_Tcell$FC_72hr_exon <- log_mean_basecsv_Tcell$Mean_72hr_exon_fpkm / log_mean_basecsv_Tcell$Mean_0hr_exon_fpkm 
#log_mean_basecsv_Tcell <- log_mean_basecsv_Tcell %>% mutate_all(funs(round(.,3)), FC_2w_exon, FC_72hr_exon)
#View(log_mean_basecsv_Tcell)

#now for the plot

ggplot(log_mean_basecsv_Tcell, aes(FC_72hr_exon,avg_ratio_72hr)) + geom_point(alpha=0.5) + ggtitle("T cell exon/intron ratio to FC: 2wks", subtitle = "applied min cutoff 2, not log transformed, replicates combined")
#NEW
ggplot(base_csv_Tcell, aes(FC_72hr_exon,Mean_72hr_ratio, label=Gene)) + geom_point(alpha=0.5) + geom_text(aes(label=ifelse(Mean_72hr_ratio>300,as.character(Gene),'')),hjust=0,vjust=0) + ggtitle("T cell exon/intron ratio to FC: 72hrs", subtitle = "applied min exon cutoff 2, not log transformed, replicated combined")

#c) ratio vs ratio fc

base_csv_Tcell$FC_72hr_ratio <- ((base_csv_Tcell$Mean_72hr_ratio - base_csv_Tcell$Mean_0hr_ratio) / base_csv_Tcell$Mean_0hr_ratio)

base_csv_Tcell$FC_72hr_ratio_LOG2 <- log2(base_csv_Tcell$FC_72hr_ratio)
base_csv_Tcell$Mean_72hr_ratio_LOG2 <- log2(base_csv_Tcell$Mean_72hr_ratio)

ggplot(base_csv_Tcell, aes(FC_72hr_ratio_LOG2,Mean_72hr_ratio_LOG2, label=Gene)) + geom_point(alpha=0.5) + geom_text_repel(aes(label=ifelse(Mean_72hr_ratio_LOG2>9,as.character(Gene),'')),hjust=0,vjust=0) + ggtitle("T cell exon/intron ratio to ratio FC: 72hrs", subtitle = "applied min exon cutoff 2, log transformed, replicates combined")


#Then more DESEQ vis stuff
      
