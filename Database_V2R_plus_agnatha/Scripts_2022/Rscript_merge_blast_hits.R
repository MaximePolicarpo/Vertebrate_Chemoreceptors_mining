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



#Now group Besthit that are less than 20000 distant

#Extend each regions 10000 on three prime and 10000 on five prime.




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


#extend all hits by 10000bp upstream and downstream

blast_rslt_col_filter_extend <- blast_rslt_col_filter %>% mutate(extend_start = start - 10000) %>% mutate(extend_end = end + 10000)
blast_rslt_col_filter_extend <- blast_rslt_col_filter_extend %>% mutate(fixed_start = case_when(
	extend_start <= 0 ~ 1,
	extend_start > 0 ~ extend_start))


blast_rslt_col_filter_extend <- blast_rslt_col_filter_extend %>% dplyr::select(seqnames, fixed_start, extend_end)
colnames(blast_rslt_col_filter_extend) <- c("seqnames", "start", "end")



#Put the table as a grange object
blast_rslt_irange <- blast_rslt_col_filter_extend %>% as_granges()


#Reduce the table to merge overlapping results
blast_rslt_disjoin <- reduce(blast_rslt_irange,with.revmap=TRUE)


#transform the table in a data.frame

Best_hits_filtered <- as.data.frame(blast_rslt_disjoin)


#Only keep regions that have a length >=100 bp
Best_hits_filtered <- Best_hits_filtered %>% filter(width >= 100)


#Merge the column with scaffold name and coordinate for samtools
Best_hits_filtered$samtools_name <- paste(Best_hits_filtered$seqnames, ":", Best_hits_filtered$start, "-", Best_hits_filtered$end, sep = "")


#Write regions in a text file
write(Best_hits_filtered$samtools_name, file="Potential_V2R_regions.tsv")







###TRASH



#list_revmap <- as.data.frame(mcols(blast_rslt_disjoin))



#filtered_data <- c()
#for(i in 1:nrow(list_revmap)){
#  filtered_data <- c(filtered_data, (slice(list_revmap, i) %>% unlist(use.names=FALSE))[which.max(slice(blast_rslt,  slice(list_revmap, i) %>% unlist(use.names=FALSE))$length)])
#}


#Best_hits_filtered <- slice(blast_rslt, filtered_data)


#write.table(potential_multiple_exon_regions, file="Potential_multiple_exon_regions.tsv", quote=FALSE, sep='\t', row.names = FALSE, col.names=FALSE)


