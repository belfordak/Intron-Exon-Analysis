#Converting mm9 gtf files to mm10

#command line work
#cat -n Exon_Only_Regions.gtf | sponge Exon_Only_Regions.gtf 
#awk -F "\t" '{print $2 "\t" $5 "\t" $6 "\t" $1 "\t" $8 "\t" $9}' Exon_Only_Regions.gtf > Exon_for_liftover.bed

#Put Exon_for_liftover.bed into UCSC liftover, output = hglft_file

#awk -F "\t" '{print $1 "," $2 "," $3 "," $4}' hglft_exon.bed > hglft_exon.csv

#awk -F "\t" '{print $1 "\t" $8 "\t" $9 "\t" $10}' Exon_Only_Regions.gtf > Exon_match_file.bed
#awk -F "\t" '{print $1 "," $2 "," $3 "," $4}' Exon_match_file.bed > Exon_match_file.csv

#import gtf as csv 
hglft_exon <- read_csv("/Users/belfordak/Data/mm9-to-mm10/hglft_exon.csv", col_names = FALSE)
match_exon <- read_csv("/Users/belfordak/Data/mm9-to-mm10/Exon_match_file.csv", col_names = FALSE)

hglft_intron <- read_csv("/Users/belfordak/Data/mm9-to-mm10/Intron/hglft_intron.csv", col_names = c("chr", "start", "end", "Num"))
match_intron <- read_csv("/Users/belfordak/Data/mm9-to-mm10/Intron/Intron_match_file.csv", col_names = c("Num", "strand", "idk", "desc"))
  
hglft_genebody <- read_csv("/Users/belfordak/Data/mm9-to-mm10/Genebody/hglft_genebody.csv", col_names = c("chr", "start", "end", "Num"))
match_genebody <- read_csv("/Users/belfordak/Data/mm9-to-mm10/Genebody/Genebody_match_file.csv", col_names = c("Num", "strand", "idk", "desc"))

#Only need to do these seperately if you say col_names=FALSE like for exons
colnames(hglft_exon) <- c("chr", "start", "end", "Num")
colnames(match_exon) <- c("Num", "strand", "idk", "desc")

#merge files by ID column
library(gdata)

merge_new_1 <- merge(hglft_exon, match_exon, by="Num") #this was just what I named my exon one
merge_intron_gtf <- merge(hglft_intron, match_intron, by="Num")
merge_genebody_gtf <- merge(hglft_genebody, match_genebody, by="Num")

#write to .csv
write.csv(merge_new_1, file = "/Users/belfordak/Data/mm9-to-mm10/Exon/exon_only_refseq_gtf.csv")
write.csv(merge_intron_gtf, file = "/Users/belfordak/Data/mm9-to-mm10/Intron/intron_only_refseq_gtf.csv")
write.csv(merge_genebody_gtf, file = "/Users/belfordak/Data/mm9-to-mm10/Genebody/genebody_refseq_gtf.csv")

#validate by namking sure $wc -l matches for this file and the hglft_.bed. They must be the same! (off by 1 though for header)
#Also, you will have two number columns from R's row.names becoming column

#next step
#txbd will be 'exonbygene'...use enxon/intron/.db in summarizeOverlap, then convert granges to gtf
