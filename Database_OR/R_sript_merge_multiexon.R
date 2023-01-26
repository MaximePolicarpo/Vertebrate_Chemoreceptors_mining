#load packages
library(data.table)
library(dplyr)
library("plyranges")
library("GenomicRanges")


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


#Extand the blast results 100bp upstream and downstream


blast_rslt <- blast_rslt %>% mutate(i_start = case_when(
  true_start > 100 ~ true_start-100,
  true_start <= 100 ~ 1
))


blast_rslt <- blast_rslt %>% mutate(i_end = true_end+100)



#Rename columns
colnames(blast_rslt) <- c("query", "seqnames", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qframe", "sframe", "qlen", "slen", "true_start", "true_end", "strand", "start", "end")


#transform the blast results table to a grange object
blast_rslt_irange <- blast_rslt %>% as_granges()

#reduce the table: merge overlapping results
blast_rslt_disjoin <- reduce(blast_rslt_irange,with.revmap=TRUE)

list_revmap <- as.data.frame(mcols(blast_rslt_disjoin))



filtered_data <- c()
for(i in 1:nrow(list_revmap)){
  filtered_data <- c(filtered_data, (slice(list_revmap, i) %>% unlist(use.names=FALSE))[which.max(slice(blast_rslt,  slice(list_revmap, i) %>% unlist(use.names=FALSE))$length)])
}


Best_hits_filtered <- slice(blast_rslt, filtered_data)


#remove regions of already found 1exon OR genes


for (row in 1:nrow(regions_functionnals_or)) {
  scaffold_name <- regions_functionnals_or[row, "scaffold"]
  debut  <-  regions_functionnals_or[row, "start"]
  fin <-  regions_functionnals_or[row, "end"]
  
  Best_hits_filtered <- Best_hits_filtered %>% filter(!(seqnames == as.character(regions_functionnals_or[row, "scaffold"]) & true_start >= (regions_functionnals_or[row, "start"] - 100) &  true_start <= (regions_functionnals_or[row, "end"] + 100)))
  Best_hits_filtered <- Best_hits_filtered %>% filter(!(seqnames == as.character(regions_functionnals_or[row, "scaffold"]) & true_end >= (regions_functionnals_or[row, "start"] - 100) &  true_end <= (regions_functionnals_or[row, "end"] + 100)))
  
}



#extand the retained regions 500bp downstream and upstream
Best_hits_filtered <- Best_hits_filtered %>% mutate(new_coord_start =  case_when(
  true_start > 501 ~ true_start-500,
  true_start < 501 ~ 1
))

Best_hits_filtered <- Best_hits_filtered %>% mutate(new_coord_end =  case_when(
  true_end < slen-501 ~ true_end+500,
  true_end > slen-501 ~ as.numeric(slen)
))


#Retain only interesting columns
Best_hits_filtered <- Best_hits_filtered %>% select(seqnames, true_start, true_end, query, new_coord_start, new_coord_end, qstart, qend)



#Now group Besthit that are less than 10000bp distant

#Extend each regions 5000bp on three prime and 5000bp on five prime. Do not extend if there is a functional genes in this windows. If so, 
#then only extend to 500bp away that the functional gene


new_coord_start_v <- c()
new_coord_end_v <- c()
for (i in seq(1:nrow(Best_hits_filtered))){

	curr_seqnames <- Best_hits_filtered[i,]$seqnames
	curr_seqnames <- as.character(curr_seqnames)
	curr_start <- Best_hits_filtered[i,]$new_coord_start
	curr_end <- Best_hits_filtered[i,]$new_coord_end
	
	if (nrow(regions_functionnals_or %>% filter(scaffold == curr_seqnames) %>% filter(end >= (curr_start-5000) & start < curr_start)) > 0){

		new_curr_start = (tail(regions_functionnals_or %>% filter(scaffold == curr_seqnames) %>% filter(end >= (curr_start-5000) & start < curr_start) %>% dplyr::arrange(end), 1)$end) + 500

	} else new_curr_start = curr_start-5000


	if (nrow(regions_functionnals_or %>% filter(scaffold == curr_seqnames) %>% filter(start <= (curr_end+5000) & end > curr_end)) > 0){

		new_curr_end = (head(regions_functionnals_or %>% filter(scaffold == curr_seqnames) %>% filter(start <= (curr_end+5000) & end > curr_end) %>% dplyr::arrange(start), 1)$start) - 500

	} else new_curr_end = curr_end+5000


	new_coord_start_v <- c(new_coord_start_v, new_curr_start)
	new_coord_end_v <- c(new_coord_end_v, new_curr_end)
	
}


Best_hits_filtered_extanded <- Best_hits_filtered %>% mutate(Extanded_start = new_coord_start_v)
Best_hits_filtered_extanded <- Best_hits_filtered_extanded %>% mutate(Extanded_end = new_coord_end_v)


# Now we use granges once again to merge overlapping best hits

colnames(Best_hits_filtered_extanded) <- c("seqnames", "true_start", "true_end", "query", "new_coord_start", "new_coord_end", "qstart", "qend", "start", "end")

Best_hits_filtered_extanded <- Best_hits_filtered_extanded %>% mutate(strand = "+")


Best_hits_filtered_irange <- Best_hits_filtered_extanded %>% as_granges()

Best_hits_filtered_disjoin <- reduce(Best_hits_filtered_irange,with.revmap=TRUE)

list_revmap <- as.data.frame(mcols(Best_hits_filtered_disjoin))




potential_multiple_exon_regions <- data.frame(NULL)
for(i in 1:nrow(list_revmap)){

	number_row <- 
	if ( length(slice(list_revmap, i) %>% unlist(use.names=FALSE)) > 0){

		line_head <- head(slice(list_revmap, i) %>% unlist(use.names=FALSE), 1)
		line_tail <- tail(slice(list_revmap, i) %>% unlist(use.names=FALSE), 1)

		start_coord <- Best_hits_filtered_extanded[line_head,]$start
		end_coord <- Best_hits_filtered_extanded[line_tail,]$end
		seqname <- as.character(Best_hits_filtered_extanded[line_head,]$seqnames)
		query <- as.character(Best_hits_filtered_extanded[line_head,]$query)

		vector_rslt <- data.frame(scaffold=seqname, start=start_coord, end=end_coord, query_id=query)

		potential_multiple_exon_regions <- rbind(potential_multiple_exon_regions, vector_rslt)
		colnames(potential_multiple_exon_regions) <- c("scaffold", "start", "end", "query_id")

	}
}


colnames(potential_multiple_exon_regions) <- c("scaffold", "start", "end", "query")


write.table(potential_multiple_exon_regions, file="Potential_multiple_exon_regions.tsv", quote=FALSE, sep='\t', row.names = FALSE, col.names=FALSE)


