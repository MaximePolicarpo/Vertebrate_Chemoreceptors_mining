#load packages

library(data.table)
library(dplyr)
library("plyranges")
library("GenomicRanges")


Best_hits_filtered <- as.data.frame(NULL) 

#Load exonerate results

exonerate_rslt <- read.table("vulgar_lines_intron_numbers_blastrslt.txt", header=FALSE, sep=" ")


#rename exonerate result columns
colnames(exonerate_rslt) <- c("query", "query_start", "query_end", "scaffold", "scaffold_start", 
                          "scaffold_end", "strand", "exonerate_score","intron_number", "best_query", "evalue")



#re-name query results for those who had no match in blastp
exonerate_rslt <- exonerate_rslt %>% mutate(best_query_F = case_when(
	best_query != "NoQuery" ~ best_query,
	best_query == "NoQuery" ~ query,
))

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
                          "i_end", "strand", "exonerate_score","intron_number", "best_query_F", "evalue")


colnames(exonerate_rslt) <- c("query", "query_start", "query_end", "scaffold", "scaffold_start", 
                          "scaffold_end", "strand", "exonerate_score","intron_number", "best_query", "evalue")



#remove exonerate results with a length <750aa. Then keep the gene that take the least amount of place. Redo the process by removing the regions of these genes


exonerate_rslt_filtered <- exonerate_rslt %>% mutate(query_length = query_end - query_start) %>% filter(query_length > 750)
exonerate_rslt_filtered <- exonerate_rslt_filtered %>% mutate(scaffold_length = scaffold_end-scaffold_start)


colnames(exonerate_rslt_filtered) <- c("query", "query_start", "query_end", "seqnames", "start", 
                          "end", "strand", "exonerate_score","intron_number", "best_query", "evalue", "query_length", "scaffold_length")



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


#Perform next step only if Best_hits_filtered is not empty


if (nrow(Best_hits_filtered) > 0) {

	#We now extend the regions 1000bp upstream and 1000bp downstream
	
	
	Best_hits_filtered <- Best_hits_filtered %>% mutate(extanded_start = start-1500) %>% mutate(extanded_end = end+1500)
	
	Best_hits_filtered <- Best_hits_filtered %>% mutate(new_coord_start =  case_when(
	  extanded_start <= 1 ~ 1,
	  extanded_start > 1 ~ extanded_start
	))
	
	
	
	#Extract the bestblast hit on these regions
	
	Best_queries_list <- c()
	for (row in 1:nrow(Best_hits_filtered)) {
	
		current_scaffold <- Best_hits_filtered[row,]$seqnames
		start_region <- Best_hits_filtered[row,]$new_coord_start
		end_region <- Best_hits_filtered[row,]$extanded_end
	
		best_query <- head(exonerate_rslt %>% mutate(query_length = query_end - query_start) %>% filter(query_length > 750) %>% filter(scaffold == current_scaffold) %>% filter(scaffold_start > start_region) %>% filter(scaffold_end < end_region) %>% arrange(evalue, desc(query_length)), 1)  %>% pull(best_query)
	
		Best_queries_list <- c(Best_queries_list, best_query)
	
	}
	
	
	Best_hits_filtered <- cbind(Best_hits_filtered, Best_queries_list)
	
	Best_hits_filtered_F <- Best_hits_filtered %>% dplyr::select(seqnames, new_coord_start, extanded_end, strand, exonerate_score, query, Best_queries_list)

} else {

	Best_hits_filtered_F <- as.data.frame(NULL) 
	Best_hits_filtered <- as.data.frame(NULL) 

}




####################################################################################
####################################################################################
####################################################################################
####################################################################################
#Now lets try to to find with a filter of 650 aa
####################################################################################
####################################################################################
####################################################################################
####################################################################################




exonerate_rslt_filtered <- exonerate_rslt %>% mutate(query_length = query_end - query_start) %>% filter(query_length > 650)
exonerate_rslt_filtered <- exonerate_rslt_filtered %>% mutate(scaffold_length = scaffold_end-scaffold_start)


colnames(exonerate_rslt_filtered) <- c("query", "query_start", "query_end", "seqnames", "start", 
                          "end", "strand", "exonerate_score","intron_number", "best_query", "evalue", "query_length", "scaffold_length")



if(nrow(Best_hits_filtered) > 0){
	for (row in 1:nrow(Best_hits_filtered)) {
		
		
		exonerate_rslt_filtered <- exonerate_rslt_filtered %>% filter(!(seqnames == as.character(Best_hits_filtered[row, "seqnames"]) & start <= (Best_hits_filtered[row, "start"]) & end >= (Best_hits_filtered[row, "start"])))
		
		
		exonerate_rslt_filtered <- exonerate_rslt_filtered %>% filter(!(seqnames == as.character(Best_hits_filtered[row, "seqnames"]) & start >= (Best_hits_filtered[row, "start"]) & end <=(Best_hits_filtered[row, "end"])))
		
		exonerate_rslt_filtered <- exonerate_rslt_filtered %>% filter(!(seqnames == as.character(Best_hits_filtered[row, "seqnames"]) & start <= (Best_hits_filtered[row, "end"]) & end >= (Best_hits_filtered[row, "end"])))
		
		exonerate_rslt_filtered <- exonerate_rslt_filtered %>% filter(!(seqnames == as.character(Best_hits_filtered[row, "seqnames"]) & start <= (Best_hits_filtered[row, "start"]) & end >= (Best_hits_filtered[row, "end"])))
		
		  
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




#We now extend the regions 1000bp upstream and 1000bp downstream



if (nrow(Best_hits_filtered_second_iter) > 0) {
	Best_hits_filtered_second_iter <- Best_hits_filtered_second_iter %>% mutate(extanded_start = start-1500) %>% mutate(extanded_end = end+1500)


	Best_hits_filtered_second_iter <- Best_hits_filtered_second_iter %>% mutate(new_coord_start =  case_when(
	  extanded_start <= 1 ~ 1,
	  extanded_start > 1 ~ extanded_start
	))



	#Extract the bestblast hit on these regions
	
	Best_queries_list <- c()
	for (row in 1:nrow(Best_hits_filtered_second_iter)) {
	
		current_scaffold <- Best_hits_filtered_second_iter[row,]$seqnames
		start_region <- Best_hits_filtered_second_iter[row,]$new_coord_start
		end_region <- Best_hits_filtered_second_iter[row,]$extanded_end
	
		best_query <- head(exonerate_rslt %>% mutate(query_length = query_end - query_start) %>% filter(query_length > 650) %>% filter(scaffold == current_scaffold) %>% filter(scaffold_start > start_region) %>% filter(scaffold_end < end_region) %>% arrange(evalue, desc(query_length)), 1)  %>% pull(best_query)
		
		Best_queries_list <- c(Best_queries_list, best_query)
	
	}
	
	
	Best_hits_filtered_second_iter <- cbind(Best_hits_filtered_second_iter, Best_queries_list)

	Best_hits_filtered_second_iter_F <- Best_hits_filtered_second_iter %>% dplyr::select(seqnames, new_coord_start, extanded_end, strand, exonerate_score, query, Best_queries_list)

} else {

	Best_hits_filtered_second_iter <- as.data.frame(NULL)
	Best_hits_filtered_second_iter_F <- as.data.frame(NULL)
}



Best_hits_filtered <- rbind(Best_hits_filtered, Best_hits_filtered_second_iter)
Best_hits_filtered_F <- rbind(Best_hits_filtered_F, Best_hits_filtered_second_iter_F)





####################################################################################
####################################################################################
####################################################################################
####################################################################################
#Now lets try to find if there are smaller genes divided in multiple exons.
####################################################################################
####################################################################################
####################################################################################
####################################################################################



exonerate_rslt_filtered <- exonerate_rslt %>% mutate(query_length = query_end - query_start) %>% filter(query_length > 200)
exonerate_rslt_filtered <- exonerate_rslt_filtered %>% mutate(scaffold_length = scaffold_end-scaffold_start)


colnames(exonerate_rslt_filtered) <- c("query", "query_start", "query_end", "seqnames", "start", 
                          "end", "strand", "exonerate_score","intron_number", "best_query", "evalue", "query_length", "scaffold_length")



if(nrow(Best_hits_filtered) > 0){
	for (row in 1:nrow(Best_hits_filtered)) {
		
		
		exonerate_rslt_filtered <- exonerate_rslt_filtered %>% filter(!(seqnames == as.character(Best_hits_filtered[row, "seqnames"]) & start <= (Best_hits_filtered[row, "start"]) & end >= (Best_hits_filtered[row, "start"])))
		
		
		exonerate_rslt_filtered <- exonerate_rslt_filtered %>% filter(!(seqnames == as.character(Best_hits_filtered[row, "seqnames"]) & start >= (Best_hits_filtered[row, "start"]) & end <=(Best_hits_filtered[row, "end"])))
		
		exonerate_rslt_filtered <- exonerate_rslt_filtered %>% filter(!(seqnames == as.character(Best_hits_filtered[row, "seqnames"]) & start <= (Best_hits_filtered[row, "end"]) & end >= (Best_hits_filtered[row, "end"])))
		
		exonerate_rslt_filtered <- exonerate_rslt_filtered %>% filter(!(seqnames == as.character(Best_hits_filtered[row, "seqnames"]) & start <= (Best_hits_filtered[row, "start"]) & end >= (Best_hits_filtered[row, "end"])))
		
		  
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




#We now extend the regions 1000bp upstream and 1000bp downstream



if (nrow(Best_hits_filtered_second_iter) > 0) {
	Best_hits_filtered_second_iter <- Best_hits_filtered_second_iter %>% mutate(extanded_start = start-1500) %>% mutate(extanded_end = end+1500)


	Best_hits_filtered_second_iter <- Best_hits_filtered_second_iter %>% mutate(new_coord_start =  case_when(
	  extanded_start <= 1 ~ 1,
	  extanded_start > 1 ~ extanded_start
	))



	#Extract the bestblast hit on these regions
	
	Best_queries_list <- c()
	for (row in 1:nrow(Best_hits_filtered_second_iter)) {
	
		current_scaffold <- Best_hits_filtered_second_iter[row,]$seqnames
		start_region <- Best_hits_filtered_second_iter[row,]$new_coord_start
		end_region <- Best_hits_filtered_second_iter[row,]$extanded_end
	
		best_query <- head(exonerate_rslt %>% mutate(query_length = query_end - query_start) %>% filter(query_length > 200) %>% filter(scaffold == current_scaffold) %>% filter(scaffold_start > start_region) %>% filter(scaffold_end < end_region) %>% arrange(evalue, desc(query_length)), 1)  %>% pull(best_query)
		
		Best_queries_list <- c(Best_queries_list, best_query)
	
	}
	
	
	Best_hits_filtered_second_iter <- cbind(Best_hits_filtered_second_iter, Best_queries_list)

	Best_hits_filtered_second_iter <- Best_hits_filtered_second_iter %>% dplyr::select(seqnames, new_coord_start, extanded_end, strand, exonerate_score, query, Best_queries_list)

} else {

	Best_hits_filtered_second_iter <- as.data.frame(NULL) 
}






#Merge the two dataframes


Best_hits_filtered <- rbind(Best_hits_filtered_F, Best_hits_filtered_second_iter)




###################################
###################################
###################################
### NEW for pseudogenes >100 <200aa


exonerate_rslt_filtered <- exonerate_rslt %>% mutate(query_length = query_end - query_start) %>% filter(query_length >= 100)
exonerate_rslt_filtered <- exonerate_rslt_filtered %>% mutate(scaffold_length = scaffold_end-scaffold_start)


colnames(exonerate_rslt_filtered) <- c("query", "query_start", "query_end", "seqnames", "start", 
                          "end", "strand", "exonerate_score","intron_number", "best_query", "evalue", "query_length", "scaffold_length")



if(nrow(Best_hits_filtered) > 0){
	for (row in 1:nrow(Best_hits_filtered)) {
		
		
		exonerate_rslt_filtered <- exonerate_rslt_filtered %>% filter(!(seqnames == as.character(Best_hits_filtered[row, "seqnames"]) & start <= (Best_hits_filtered[row, "new_coord_start"]) & end >= (Best_hits_filtered[row, "new_coord_start"])))
		
		
		exonerate_rslt_filtered <- exonerate_rslt_filtered %>% filter(!(seqnames == as.character(Best_hits_filtered[row, "seqnames"]) & start >= (Best_hits_filtered[row, "new_coord_start"]) & end <=(Best_hits_filtered[row, "extanded_end"])))
		
		exonerate_rslt_filtered <- exonerate_rslt_filtered %>% filter(!(seqnames == as.character(Best_hits_filtered[row, "seqnames"]) & start <= (Best_hits_filtered[row, "extanded_end"]) & end >= (Best_hits_filtered[row, "extanded_end"])))
		
		exonerate_rslt_filtered <- exonerate_rslt_filtered %>% filter(!(seqnames == as.character(Best_hits_filtered[row, "seqnames"]) & start <= (Best_hits_filtered[row, "new_coord_start"]) & end >= (Best_hits_filtered[row, "extanded_end"])))
		
		  
	}
}




Best_hits_filtered_third_iter <- as.data.frame(NULL) 



while (nrow(exonerate_rslt_filtered) > 0){


	exonerate_rslt_filtered_irange <- exonerate_rslt_filtered %>% as_granges()
	exonerate_rslt_filtered_disjoin <- reduce(exonerate_rslt_filtered_irange,with.revmap=TRUE)
	list_revmap <- as.data.frame(mcols(exonerate_rslt_filtered_disjoin))

	filtered_data <- c()
	for(i in 1:nrow(list_revmap)){
  	filtered_data <- c(filtered_data, (slice(list_revmap, i) %>% unlist(use.names=FALSE))[which.min(slice(exonerate_rslt_filtered, slice(list_revmap, i) %>% unlist(use.names=FALSE))$scaffold_length)])
	}


	Best_hits_filtered_third_iter <- rbind(Best_hits_filtered_third_iter, slice(exonerate_rslt_filtered, filtered_data))


	if(nrow(Best_hits_filtered_third_iter) > 0){
		for (row in 1:nrow(Best_hits_filtered_third_iter)) {
		
		
			exonerate_rslt_filtered <- exonerate_rslt_filtered %>% filter(!(seqnames == as.character(Best_hits_filtered_third_iter[row, "seqnames"]) & start <= (	Best_hits_filtered_third_iter[row, "start"]) & end >= (Best_hits_filtered_third_iter[row, "start"])))
		
		
			exonerate_rslt_filtered <- exonerate_rslt_filtered %>% filter(!(seqnames == as.character(Best_hits_filtered_third_iter[row, "seqnames"]) & start >= (	Best_hits_filtered_third_iter[row, "start"]) & end <=(Best_hits_filtered_third_iter[row, "end"])))
		
			exonerate_rslt_filtered <- exonerate_rslt_filtered %>% filter(!(seqnames == as.character(Best_hits_filtered_third_iter[row, "seqnames"]) & start <= (	Best_hits_filtered_third_iter[row, "end"]) & end >= (Best_hits_filtered_third_iter[row, "end"])))
		
			exonerate_rslt_filtered <- exonerate_rslt_filtered %>% filter(!(seqnames == as.character(Best_hits_filtered_third_iter[row, "seqnames"]) & start <= (	Best_hits_filtered_third_iter[row, "start"]) & end >= (Best_hits_filtered_third_iter[row, "end"])))
		
		  
		}
	
	}



}




#We now extend the regions 100bp upstream and 100bp downstream


if (nrow(Best_hits_filtered_third_iter) > 0) {

	Best_hits_filtered_third_iter <- Best_hits_filtered_third_iter %>% mutate(extanded_start = start-100) %>% mutate(extanded_end = end+100)
	
	
	Best_hits_filtered_third_iter <- Best_hits_filtered_third_iter %>% mutate(new_coord_start =  case_when(
	  extanded_start <= 1 ~ 1,
	  extanded_start > 1 ~ extanded_start
	))
	
	
	
	#Extract the bestblast hit on these regions
	
	Best_queries_list <- c()
	for (row in 1:nrow(Best_hits_filtered_third_iter)) {
	
		current_scaffold <- Best_hits_filtered_third_iter[row,]$seqnames
		start_region <- Best_hits_filtered_third_iter[row,]$new_coord_start
		end_region <- Best_hits_filtered_third_iter[row,]$extanded_end
	
		best_query <- head(exonerate_rslt %>% mutate(query_length = query_end - query_start) %>% filter(query_length >= 100) %>% filter(scaffold == current_scaffold) %>% filter(scaffold_start > start_region) %>% filter(scaffold_end < end_region) %>% arrange(evalue, desc(query_length)), 1) %>% pull(best_query)
	
		Best_queries_list <- c(Best_queries_list, best_query)
	
	}
	
	
	Best_hits_filtered_third_iter <- cbind(Best_hits_filtered_third_iter, Best_queries_list)
	
	Best_hits_filtered_third_iter <- Best_hits_filtered_third_iter %>% dplyr::select(seqnames, new_coord_start, extanded_end, strand, exonerate_score, query, Best_queries_list)

} else {

	Best_hits_filtered_third_iter <- as.data.frame(NULL) 
}




#write the result in a table
write.table(Best_hits_filtered_third_iter, file="Parsed_exonerate_pseudogene_regions.tsv", quote=FALSE, sep='\t', row.names = FALSE, col.names=FALSE)








