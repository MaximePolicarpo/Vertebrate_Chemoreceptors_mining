library(data.table)
library(dplyr)
library("plyranges")
library("GenomicRanges")

blast_rslt <- read.table("tblastn_functionnal_or_vs_genome.tblastn", header=FALSE, sep="\t")
colnames(blast_rslt) <- c("query", "sseqid", "pident", "length", "mismatch", 
                          "gapopen", "qstart", "qend", "sstart", "send",
                          "evalue", "bitscore", "qframe", "sframe", "qlen", "slen")



#Extract regions of functionnals olfactory receptors


regions_functionnals_or <- read.table("Coordinates_Functionnal_ORS.txt", header=FALSE, sep="\t")

colnames(regions_functionnals_or) <- c("scaffold", "start", "end")


blast_rslt <- blast_rslt %>% mutate(true_start = case_when(
  sstart < send ~ sstart,
  send < sstart ~ send
))


blast_rslt <- blast_rslt %>% mutate(true_end = case_when(
  sstart < send ~ send,
  send < sstart ~ sstart
))


blast_rslt <- blast_rslt %>% mutate(strand = case_when(
  sstart < send ~ "+",
  send < sstart ~ "-"
))


##


blast_rslt <- blast_rslt %>% mutate(i_start = case_when(
  true_start > 100 ~ true_start-100,
  true_start <= 100 ~ as.numeric(true_start)
))


blast_rslt <- blast_rslt %>% mutate(i_end = true_end+100)




##





colnames(blast_rslt) <- c("query", "seqnames", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qframe", "sframe", "qlen", "slen", "true_start", "true_end", "strand", "start", "end")



blast_rslt_irange <- blast_rslt %>% as_granges()

blast_rslt_disjoin <- reduce(blast_rslt_irange,with.revmap=TRUE)

list_revmap <- as.data.frame(mcols(blast_rslt_disjoin))



filtered_data <- c()
for(i in 1:nrow(list_revmap)){
  filtered_data <- c(filtered_data, (slice(list_revmap, i) %>% unlist(use.names=FALSE))[which.max(slice(blast_rslt,  slice(list_revmap, i) %>% unlist(use.names=FALSE))$length)])
}


Best_hits_filtered <- slice(blast_rslt, filtered_data)


#remove regions of functionnals#


for (row in 1:nrow(regions_functionnals_or)) {
  scaffold_name <- regions_functionnals_or[row, "scaffold"]
  debut  <-  regions_functionnals_or[row, "start"]
  fin <-  regions_functionnals_or[row, "end"]
  
  Best_hits_filtered <- Best_hits_filtered %>% filter(!(seqnames == as.character(regions_functionnals_or[row, "scaffold"]) & true_start >= (regions_functionnals_or[row, "start"] - 100) &  true_start <= (regions_functionnals_or[row, "end"] + 100)))
  Best_hits_filtered <- Best_hits_filtered %>% filter(!(seqnames == as.character(regions_functionnals_or[row, "scaffold"]) & true_end >= (regions_functionnals_or[row, "start"] - 100) &  true_end <= (regions_functionnals_or[row, "end"] + 100)))
  
}


Best_hits_filtered <- Best_hits_filtered %>% mutate(new_coord_start =  case_when(
  true_start > 501 ~ true_start-500,
  true_start < 501 ~ 1
))

Best_hits_filtered <- Best_hits_filtered %>% mutate(new_coord_end =  case_when(
  true_end < slen-501 ~ true_end+500,
  true_end > slen-501 ~ as.numeric(slen)
))


Best_hits_filtered <- Best_hits_filtered %>% select(seqnames, true_start, true_end, query, new_coord_start, new_coord_end, qstart, qend)

write.table(Best_hits_filtered, file='Pseudo_truncated_coordinates.tsv', quote=FALSE, sep='\t', row.names = FALSE,
            col.names=FALSE)




