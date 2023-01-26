#load packages
library(data.table)
library(dplyr)
library("plyranges")
library("GenomicRanges")

args = commandArgs(trailingOnly=TRUE)
extension_length <- as.numeric(args[1])

#load tblastn results
blast_rslt <- read.table("tblastn_functionnal_and_known_or_vs_genome.tblastn", header=FALSE, sep="\t")

#rename column of the blast result table
colnames(blast_rslt) <- c("query", "sseqid", "pident", "length", "mismatch", 
                          "gapopen", "qstart", "qend", "sstart", "send",
                          "evalue", "bitscore", "qframe", "sframe", "qlen", "slen")



#Extract regions of already found 1exon functionnals olfactory receptors
regions_functionnals_or <- read.table("Coordinates_Functionnal_ORS.txt", header=FALSE, sep="\t")

#name columns
colnames(regions_functionnals_or) <- c("scaffold", "start", "end")


#re-orientate blast results
blast_rslt <- blast_rslt %>% mutate(true_start = case_when(
  sstart < send ~ sstart,
  send < sstart ~ send
))


blast_rslt <- blast_rslt %>% mutate(true_end = case_when(
  sstart < send ~ send,
  send < sstart ~ sstart
))


#Tell in which strand results are
blast_rslt <- blast_rslt %>% mutate(strand = case_when(
  sstart < send ~ "+",
  send < sstart ~ "-"
))


#select columns
blast_rslt <- blast_rslt %>% dplyr::select("query", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "evalue", "bitscore", "qframe", "sframe", "qlen", "slen", "true_start", "true_end", "strand")


#Rename columns
colnames(blast_rslt) <- c("query", "seqnames", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "evalue", "bitscore", "qframe", "sframe", "qlen", "slen", "start", "end", "strand")


#transform the blast results table to a grange object
blast_rslt_irange <- blast_rslt %>% as_granges()

#reduce the table: merge overlapping results
blast_rslt_disjoin <- reduce(blast_rslt_irange,with.revmap=TRUE)
Best_hits_filtered <- as.data.frame(blast_rslt_disjoin)


#remove regions of already found 1exon OR genes


for (row in 1:nrow(regions_functionnals_or)) {

  Best_hits_filtered <- Best_hits_filtered %>% filter(!(seqnames == as.character(regions_functionnals_or[row, "scaffold"]) & start <= (regions_functionnals_or[row, "start"]) & end >= (regions_functionnals_or[row, "start"])))
  Best_hits_filtered <- Best_hits_filtered %>% filter(!(seqnames == as.character(regions_functionnals_or[row, "scaffold"]) & start >= (regions_functionnals_or[row, "start"]) & end <= (regions_functionnals_or[row, "end"])))
  Best_hits_filtered <- Best_hits_filtered %>% filter(!(seqnames == as.character(regions_functionnals_or[row, "scaffold"]) & start <= (regions_functionnals_or[row, "end"]) & end >= (regions_functionnals_or[row, "end"])))


}



#Now group Besthit that are less than 10000bp distant

#Extend each regions 5000bp on three prime and 5000bp on five prime. Do not extend if there is a functional genes in this windows. If so, 
#then only extend to 50bp away that the functional gene


new_coord_start_v <- c()
new_coord_end_v <- c()
for (i in seq(1:nrow(Best_hits_filtered))){

	curr_seqnames <- Best_hits_filtered[i,]$seqnames
	curr_seqnames <- as.character(curr_seqnames)
	curr_start <- Best_hits_filtered[i,]$start
	curr_end <- Best_hits_filtered[i,]$end
	
	if (nrow(regions_functionnals_or %>% filter(scaffold == curr_seqnames) %>% filter(end >= (curr_start-extension_length) & start < curr_start)) > 0){

		new_curr_start = (tail(regions_functionnals_or %>% filter(scaffold == curr_seqnames) %>% filter(end >= (curr_start-extension_length) & start < curr_start) %>% dplyr::arrange(end), 1)$end) + 50

	} else new_curr_start = curr_start-extension_length


	if (nrow(regions_functionnals_or %>% filter(scaffold == curr_seqnames) %>% filter(start <= (curr_end+extension_length) & end > curr_end)) > 0){

		new_curr_end = (head(regions_functionnals_or %>% filter(scaffold == curr_seqnames) %>% filter(start <= (curr_end+extension_length) & end > curr_end) %>% dplyr::arrange(start), 1)$start) - 50

	} else new_curr_end = curr_end+extension_length


	new_coord_start_v <- c(new_coord_start_v, new_curr_start)
	new_coord_end_v <- c(new_coord_end_v, new_curr_end)
	
}


Best_hits_filtered_extanded <- Best_hits_filtered %>% mutate(Extanded_start = new_coord_start_v)
Best_hits_filtered_extanded <- Best_hits_filtered_extanded %>% mutate(Extanded_end = new_coord_end_v)




#Retain only interesting columns and put 1 if the start coordinate is negative
Best_hits_filtered_extanded <- Best_hits_filtered_extanded %>% mutate(Extanded_start_fixed =  case_when(
  Extanded_start <= 0 ~ 1,
  Extanded_start > 0 ~ Extanded_start
))



Best_hits_filtered_extanded <- Best_hits_filtered_extanded %>% select(seqnames, Extanded_start_fixed, Extanded_end)
colnames(Best_hits_filtered_extanded) <- c("seqnames", "start", "end")


# Now we use granges once again to merge overlapping extanded best hits

Best_hits_filtered_irange <- Best_hits_filtered_extanded %>% as_granges()
Best_hits_filtered_disjoin <- reduce(Best_hits_filtered_irange,with.revmap=TRUE)
Best_hits_filtered_extanded_joined <- as.data.frame(Best_hits_filtered_disjoin)


#Merge the column with scaffold name and coordinate for samtools
Best_hits_filtered_extanded_joined$samtools_name <- paste(Best_hits_filtered_extanded_joined$seqnames, ":", Best_hits_filtered_extanded_joined$start, "-", Best_hits_filtered_extanded_joined$end, sep = "")


#Write regions in a text file
write(Best_hits_filtered_extanded_joined$samtools_name, file="Potential_multiple_exon_regions.tsv")







###TRASH



#list_revmap <- as.data.frame(mcols(blast_rslt_disjoin))



#filtered_data <- c()
#for(i in 1:nrow(list_revmap)){
#  filtered_data <- c(filtered_data, (slice(list_revmap, i) %>% unlist(use.names=FALSE))[which.max(slice(blast_rslt,  slice(list_revmap, i) %>% unlist(use.names=FALSE))$length)])
#}


#Best_hits_filtered <- slice(blast_rslt, filtered_data)


#write.table(potential_multiple_exon_regions, file="Potential_multiple_exon_regions.tsv", quote=FALSE, sep='\t', row.names = FALSE, col.names=FALSE)


