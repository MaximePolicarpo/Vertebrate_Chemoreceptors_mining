#load packages
library(data.table)
library(dplyr)
library("plyranges")
library("GenomicRanges")


#load tblastn results
blast_rslt <- read.table("V2R_vs_Genome.blastn", header=FALSE, sep="\t")

#rename column of the blast result table

colnames(blast_rslt) <- c("query", "sseqid", "pident", "length", "mismatch", 
                          "gapopen", "qstart", "qend", "sstart", "send",
                          "evalue", "bitscore")



#extract interesting columns
blast_rslt_col_filter <- blast_rslt %>% dplyr::select(sseqid, sstart, send, evalue, length) 


#re-orientate results
blast_rslt_col_filter <- blast_rslt_col_filter %>% mutate(i_start = case_when(
  sstart < send ~ sstart,
  send < sstart ~ send
))


blast_rslt_col_filter <- blast_rslt_col_filter %>% mutate(i_end = case_when(
  sstart < send ~ send,
  send < sstart ~ sstart
))

#put a column to indicate on which strand the gene is
blast_rslt_col_filter <- blast_rslt_col_filter %>% mutate(strand = case_when(
  sstart < send ~ "+",
  send < sstart ~ "-"
))



#We select only interesting columns with the new start and end re-orientated
blast_rslt_col_filter <- blast_rslt_col_filter %>% select(sseqid, i_start, i_end, evalue, strand, length)

#rename columns
colnames(blast_rslt_col_filter) <- c("seqnames", "start", "end", "evalue", "strand", "length")


#Put the table as a grange object
blast_rslt_irange <- blast_rslt_col_filter %>% as_granges()


#Reduce the table to merge overlapping results
blast_rslt_disjoin <- reduce(blast_rslt_irange,with.revmap=TRUE)


#transform the table in a data.frame

Best_hits_filtered <- as.data.frame(blast_rslt_disjoin)


#Only keep regions that have a length >=100 bp
Best_hits_filtered <- Best_hits_filtered %>% filter(width >= 50)


#Merge the column with scaffold name and coordinate for samtools
Best_hits_filtered$samtools_name <- paste(Best_hits_filtered$seqnames, ":", Best_hits_filtered$start, "-", Best_hits_filtered$end, sep = "")


V2R_best_hits <- scan(file="V2R_Regions.tsv", what="character")


Best_hits_filtered <- Best_hits_filtered %>% filter(samtools_name %in% V2R_best_hits)



#extend all hits by 10000bp upstream and downstream

Best_hits_filtered_extend <- Best_hits_filtered %>% mutate(extend_start = start - 10000) %>% mutate(extend_end = end + 10000)
Best_hits_filtered_extend <- Best_hits_filtered_extend %>% mutate(fixed_start = case_when(
	extend_start <= 0 ~ 1,
	extend_start > 0 ~ extend_start))


Best_hits_filtered_extend <- Best_hits_filtered_extend %>% dplyr::select(seqnames, fixed_start, extend_end)
colnames(Best_hits_filtered_extend) <- c("seqnames", "start", "end")



#Put the table as a grange object
Best_hits_irange <- Best_hits_filtered_extend %>% as_granges()


#Reduce the table to merge overlapping results
Best_hits_disjoin <- reduce(Best_hits_irange,with.revmap=TRUE)


#transform the table in a data.frame

Best_hits_filtered <- as.data.frame(Best_hits_disjoin)


#Only keep regions that have a length >=100 bp
Best_hits_filtered <- Best_hits_filtered %>% filter(width >= 100)


#Merge the column with scaffold name and coordinate for samtools
Best_hits_filtered$samtools_name <- paste(Best_hits_filtered$seqnames, ":", Best_hits_filtered$start, "-", Best_hits_filtered$end, sep = "")


#Write regions in a text file

write(Best_hits_filtered$samtools_name, file="Potential_V2R_regions.tsv")













