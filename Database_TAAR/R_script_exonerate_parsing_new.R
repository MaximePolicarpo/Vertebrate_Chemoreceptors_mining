#load packages

library(data.table)
library(dplyr)
library("plyranges")
library("GenomicRanges")


Best_hits_filtered <- as.data.frame(NULL) 

#Load exonerate results

exonerate_rslt <- read.table("vulgar_lines_intron_numbers.txt", header=FALSE, sep=" ")


#rename exonerate result columns
colnames(exonerate_rslt) <- c("query", "query_start", "query_end", "scaffold", "scaffold_start", 
                          "scaffold_end", "strand", "exonerate_score","intron_number")



#re-orientate exonerate results coordinates


exonerate_rslt <- exonerate_rslt %>% mutate(i_start = case_when(
  scaffold_start < scaffold_end ~ scaffold_start,
  scaffold_end < scaffold_start ~ scaffold_end
))


exonerate_rslt <- exonerate_rslt %>% mutate(i_end = case_when(
  scaffold_start < scaffold_end ~ scaffold_end,
  scaffold_end < scaffold_start ~ scaffold_start
))



exonerate_rslt <- exonerate_rslt %>% dplyr::select("query", "query_start", "query_end", "scaffold", "i_start", 
                          "i_end", "strand", "exonerate_score","intron_number")


colnames(exonerate_rslt) <- c("query", "query_start", "query_end", "scaffold", "scaffold_start", 
                          "scaffold_end", "strand", "exonerate_score","intron_number")








#remove exonerate results with less than one intron and <230aa. Then keep the gene that take the least amount of place. Redo the process by removing the regions of these genes


exonerate_rslt_filtered <- exonerate_rslt %>% mutate(query_length = query_end - query_start) %>% filter(intron_number >= 1) %>% filter(query_length > 250)
exonerate_rslt_filtered <- exonerate_rslt_filtered %>% mutate(scaffold_length = scaffold_end-scaffold_start)


colnames(exonerate_rslt_filtered) <- c("query", "query_start", "query_end", "seqnames", "start", 
                          "end", "strand", "exonerate_score","intron_number", "query_length", "scaffold_length")





while (nrow(exonerate_rslt_filtered) > 0){


	exonerate_rslt_filtered_irange <- exonerate_rslt_filtered %>% as_granges()
	exonerate_rslt_filtered_disjoin <- reduce(exonerate_rslt_filtered_irange,with.revmap=TRUE)
	list_revmap <- as.data.frame(mcols(exonerate_rslt_filtered_disjoin))

	filtered_data <- c()
	for(i in 1:nrow(list_revmap)){
  	filtered_data <- c(filtered_data, (slice(list_revmap, i) %>% unlist(use.names=FALSE))[which.min(slice(exonerate_rslt_filtered, slice(list_revmap, i) %>% unlist(use.names=FALSE))$scaffold_length)])
	}


	Best_hits_filtered <- rbind(Best_hits_filtered, slice(exonerate_rslt_filtered, filtered_data))


	if(nrow(Best_hits_filtered) > 0){
		for (row in 1:nrow(Best_hits_filtered)) {
		
		
			exonerate_rslt_filtered <- exonerate_rslt_filtered %>% filter(!(seqnames == as.character(Best_hits_filtered[row, "seqnames"]) & start <= (	Best_hits_filtered[row, "start"]) & end >= (Best_hits_filtered[row, "start"])))
		
		
			exonerate_rslt_filtered <- exonerate_rslt_filtered %>% filter(!(seqnames == as.character(Best_hits_filtered[row, "seqnames"]) & start >= (	Best_hits_filtered[row, "start"]) & end <=(Best_hits_filtered[row, "end"])))
		
			exonerate_rslt_filtered <- exonerate_rslt_filtered %>% filter(!(seqnames == as.character(Best_hits_filtered[row, "seqnames"]) & start <= (	Best_hits_filtered[row, "end"]) & end >= (Best_hits_filtered[row, "end"])))
		
			exonerate_rslt_filtered <- exonerate_rslt_filtered %>% filter(!(seqnames == as.character(Best_hits_filtered[row, "seqnames"]) & start <= (	Best_hits_filtered[row, "start"]) & end >= (Best_hits_filtered[row, "end"])))
		
		  
		}
	
	}



}


#We now extend the regions 800bp upstream and 800bp downstream

Best_hits_filtered <- Best_hits_filtered %>% mutate(extanded_start = start-800) %>% mutate(extanded_end = end+800)

Best_hits_filtered <- Best_hits_filtered %>% mutate(new_coord_start =  case_when(
  extanded_start <= 1 ~ 1,
  extanded_start > 1 ~ extanded_start
))




Best_hits_filtered_F <- Best_hits_filtered %>% dplyr::select(seqnames, new_coord_start, extanded_end, strand, exonerate_score, query)


#write the result in a table
write.table(Best_hits_filtered_F, file="Filtered_Potential_multiple_exon_regions.tsv", quote=FALSE, sep='\t', row.names = FALSE, col.names=FALSE)

######
######
######
######
######
######
#Now lets try to find if there are smaller genes divided in multiple exons.




exonerate_rslt_filtered <- exonerate_rslt %>% mutate(query_length = query_end - query_start) %>% filter(intron_number >= 1) %>% filter(query_length > 200)
exonerate_rslt_filtered <- exonerate_rslt_filtered %>% mutate(scaffold_length = scaffold_end-scaffold_start)


colnames(exonerate_rslt_filtered) <- c("query", "query_start", "query_end", "seqnames", "start", 
                          "end", "strand", "exonerate_score","intron_number", "query_length", "scaffold_length")



if(nrow(Best_hits_filtered) > 0){
	for (row in 1:nrow(Best_hits_filtered)) {
		
		
		exonerate_rslt_filtered <- exonerate_rslt_filtered %>% filter(!(seqnames == as.character(Best_hits_filtered[row, "seqnames"]) & start <= (	Best_hits_filtered[row, "start"]) & end >= (Best_hits_filtered[row, "start"])))
		
		
		exonerate_rslt_filtered <- exonerate_rslt_filtered %>% filter(!(seqnames == as.character(Best_hits_filtered[row, "seqnames"]) & start >= (	Best_hits_filtered[row, "start"]) & end <=(Best_hits_filtered[row, "end"])))
		
		exonerate_rslt_filtered <- exonerate_rslt_filtered %>% filter(!(seqnames == as.character(Best_hits_filtered[row, "seqnames"]) & start <= (	Best_hits_filtered[row, "end"]) & end >= (Best_hits_filtered[row, "end"])))
		
		exonerate_rslt_filtered <- exonerate_rslt_filtered %>% filter(!(seqnames == as.character(Best_hits_filtered[row, "seqnames"]) & start <= (	Best_hits_filtered[row, "start"]) & end >= (Best_hits_filtered[row, "end"])))
		
		  
	}
}




Best_hits_filtered_second_iter <- as.data.frame(NULL) 



while (nrow(exonerate_rslt_filtered) > 0){


	exonerate_rslt_filtered_irange <- exonerate_rslt_filtered %>% as_granges()
	exonerate_rslt_filtered_disjoin <- reduce(exonerate_rslt_filtered_irange,with.revmap=TRUE)
	list_revmap <- as.data.frame(mcols(exonerate_rslt_filtered_disjoin))

	filtered_data <- c()
	for(i in 1:nrow(list_revmap)){
  	filtered_data <- c(filtered_data, (slice(list_revmap, i) %>% unlist(use.names=FALSE))[which.min(slice(exonerate_rslt_filtered, slice(list_revmap, i) %>% unlist(use.names=FALSE))$scaffold_length)])
	}


	Best_hits_filtered_second_iter <- rbind(Best_hits_filtered_second_iter, slice(exonerate_rslt_filtered, filtered_data))


	if(nrow(Best_hits_filtered_second_iter) > 0){
		for (row in 1:nrow(Best_hits_filtered_second_iter)) {
		
		
			exonerate_rslt_filtered <- exonerate_rslt_filtered %>% filter(!(seqnames == as.character(Best_hits_filtered_second_iter[row, "seqnames"]) & start <= (	Best_hits_filtered_second_iter[row, "start"]) & end >= (Best_hits_filtered_second_iter[row, "start"])))
		
		
			exonerate_rslt_filtered <- exonerate_rslt_filtered %>% filter(!(seqnames == as.character(Best_hits_filtered_second_iter[row, "seqnames"]) & start >= (	Best_hits_filtered_second_iter[row, "start"]) & end <=(Best_hits_filtered_second_iter[row, "end"])))
		
			exonerate_rslt_filtered <- exonerate_rslt_filtered %>% filter(!(seqnames == as.character(Best_hits_filtered_second_iter[row, "seqnames"]) & start <= (	Best_hits_filtered_second_iter[row, "end"]) & end >= (Best_hits_filtered_second_iter[row, "end"])))
		
			exonerate_rslt_filtered <- exonerate_rslt_filtered %>% filter(!(seqnames == as.character(Best_hits_filtered_second_iter[row, "seqnames"]) & start <= (	Best_hits_filtered_second_iter[row, "start"]) & end >= (Best_hits_filtered_second_iter[row, "end"])))
		
		  
		}
	
	}



}


#We now extend the regions 800bp upstream and 800bp downstream

Best_hits_filtered_second_iter <- Best_hits_filtered_second_iter %>% mutate(extanded_start = start-800) %>% mutate(extanded_end = end+800)


Best_hits_filtered_second_iter <- Best_hits_filtered_second_iter %>% mutate(new_coord_start =  case_when(
  extanded_start <= 1 ~ 1,
  extanded_start > 1 ~ extanded_start
))


Best_hits_filtered_second_iter <- Best_hits_filtered_second_iter %>% dplyr::select(seqnames, new_coord_start, extanded_end, strand, exonerate_score, query)


#write the result in a table
write.table(Best_hits_filtered_second_iter, file="PseudoTrunc_Potential_multiple_exon_regions.tsv", quote=FALSE, sep='\t', row.names = FALSE, col.names=FALSE)










