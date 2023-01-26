#load packages
library(data.table)
library(dplyr)
library("plyranges")
library("GenomicRanges")

args = commandArgs(trailingOnly=TRUE)
extension_length <- as.numeric(args[1])

#load tblastn results
blast_rslt <- read.table("V2R_vs_Genome.blastn", header=FALSE, sep="\t")

#rename column of the blast result table

colnames(blast_rslt) <- c("query", "sseqid", "pident", "length", "mismatch", 
                          "gapopen", "qstart", "qend", "sstart", "send",
                          "evalue", "bitscore")




#load and rename column of regions with found genes

found_genes <- read.table("Coordinates_already_examined.tsv", header=FALSE, sep="\t")

colnames(found_genes) <- c("scaffold", "start", "end")




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


#remove regions of already found V2R genes
for (row in 1:nrow(found_genes)) {

  Best_hits_filtered <- Best_hits_filtered %>% filter(!(seqnames == as.character(found_genes[row, "scaffold"]) & start <= (found_genes[row, "start"]) & end >= (found_genes[row, "start"])))
  Best_hits_filtered <- Best_hits_filtered %>% filter(!(seqnames == as.character(found_genes[row, "scaffold"]) & start >= (found_genes[row, "start"]) & end <= (found_genes[row, "end"])))
  Best_hits_filtered <- Best_hits_filtered %>% filter(!(seqnames == as.character(found_genes[row, "scaffold"]) & start <= (found_genes[row, "end"]) & end >= (found_genes[row, "end"])))


}



#Extend each regions 10000bp on three prime and 10000bp on five prime. Do not extend if there is a functional genes in this windows. If so, 
#then only extend to 50bp away that the functional gene


if (nrow(Best_hits_filtered) > 0) {
  new_coord_start_v <- c()
  new_coord_end_v <- c()
  for (i in seq(1:nrow(Best_hits_filtered))){
  
    curr_seqnames <- Best_hits_filtered[i,]$seqnames
    curr_seqnames <- as.character(curr_seqnames)
    curr_start <- Best_hits_filtered[i,]$start
    curr_end <- Best_hits_filtered[i,]$end
    
    if (nrow(found_genes %>% filter(scaffold == curr_seqnames) %>% filter(end >= (curr_start-extension_length) & start < curr_start)) > 0){
  
      new_curr_start = (tail(found_genes %>% filter(scaffold == curr_seqnames) %>% filter(end >= (curr_start-extension_length) & start < curr_start) %>% dplyr::arrange(end), 1)$end) + 50
  
    } else new_curr_start = curr_start-extension_length
  
  
    if (nrow(found_genes %>% filter(scaffold == curr_seqnames) %>% filter(start <= (curr_end+extension_length) & end > curr_end)) > 0){
  
      new_curr_end = (head(found_genes %>% filter(scaffold == curr_seqnames) %>% filter(start <= (curr_end+extension_length) & end > curr_end) %>% dplyr::arrange(start), 1)$start) - 50
  
    } else new_curr_end = curr_end+extension_length
  
  
    new_coord_start_v <- c(new_coord_start_v, new_curr_start)
    new_coord_end_v <- c(new_coord_end_v, new_curr_end)
    
  }
  
  
  Best_hits_filtered_extanded <- Best_hits_filtered %>% mutate(Extanded_start = new_coord_start_v)
  Best_hits_filtered_extanded <- Best_hits_filtered_extanded %>% mutate(Extanded_end = new_coord_end_v)
  
  
  
  
  Best_hits_filtered <- Best_hits_filtered_extanded %>% dplyr::select(seqnames, Extanded_start, Extanded_end)
  colnames(Best_hits_filtered) <- c("seqnames", "start", "end")
  
  
  
  #Put the table as a grange object
  Best_hits_filtered_irange <- Best_hits_filtered %>% as_granges()
  
  
  #Reduce the table to merge overlapping results
  Best_hits_filtered_disjoin <- reduce(Best_hits_filtered_irange,with.revmap=TRUE)
  
  
  #transform the table in a data.frame
  
  Best_hits_filtered <- as.data.frame(Best_hits_filtered_disjoin)
  
  
  #Only keep regions that have a length >=100 bp
  Best_hits_filtered <- Best_hits_filtered %>% filter(width >= 100)
  
  
  #put 1 if start coord is below 1
  Best_hits_filtered <- Best_hits_filtered %>% mutate(fixed_start = case_when(
    as.numeric(start) <= 0 ~ 1,
    as.numeric(start) > 0 ~ as.numeric(start)))
  
  
  #Merge the column with scaffold name and coordinate for samtools
  Best_hits_filtered$samtools_name <- paste(Best_hits_filtered$seqnames, ":", Best_hits_filtered$fixed_start, "-", Best_hits_filtered$end, sep = "")

} else {

	Best_hits_filtered$samtools_name <- c()
}


#Write regions in a text file
write(Best_hits_filtered$samtools_name, file="Potential_V2R_regions.tsv")






