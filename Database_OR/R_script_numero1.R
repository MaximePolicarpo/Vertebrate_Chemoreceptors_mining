#load packages

library(data.table)
library(dplyr)
library("plyranges")
library("GenomicRanges")


#Load blast result

blast_rslt <- read.table("OR_vs_Genome.blastn", header=FALSE, sep="\t")

#rename blast result columns
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
Best_hits_filtered <- Best_hits_filtered %>% filter(width >= 100)


#Extand results 1000bp upstream and 1000bp downstream

Best_hits_filtered <- Best_hits_filtered %>% mutate(new_coord_start =  case_when(
  start >= 1001 ~ start-1000,
  start < 1001 ~ 1
))

Best_hits_filtered <- Best_hits_filtered %>% mutate(new_coord_end =  end+1000)



#Merge the column with scaffold name and coordinate for samtools
Best_hits_filtered$samtools_name <- paste(Best_hits_filtered$seqnames, ":", Best_hits_filtered$new_coord_start, "-", Best_hits_filtered$new_coord_end, sep = "")


#Write regions in a text file
write(Best_hits_filtered$samtools_name, file="Best_hits_filtered.tsv")




#########
#########
##TRASH##
#########
#########
#########

#list_revmap <- as.data.frame(mcols(blast_rslt_disjoin))



#filtered_data <- c()
#for(i in 1:nrow(list_revmap)){
#  filtered_data <- c(filtered_data, (slice(list_revmap, i) %>% unlist(use.names=FALSE))[which.min(slice(blast_rslt_col_filter,  slice(list_revmap, i) %>% unlist(use.names=FALSE))$evalue)])
#}


#Best_hits_filtered <- slice(blast_rslt_col_filter, filtered_data)
#Best_hits_filtered <- Best_hits_filtered %>% filter(length >= 100)




#Keep only interesting columsn to write in the final result file
#Best_hits_filtered_F <- Best_hits_filtered %>% select(seqnames, new_coord_start, new_coord_end)


#write the final result file
#write.table(Best_hits_filtered_F, file='Best_hits_filtered.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names=FALSE)




