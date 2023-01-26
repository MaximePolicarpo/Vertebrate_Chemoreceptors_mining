#load packages

library(data.table)
library(dplyr)
library("plyranges")
library("GenomicRanges")
library("tidyr")
library(stringr)



#Load exonerate results

exonerate_rslt <- read.table("vulgar_lines_intron_numbers_blastrslt.txt", header=FALSE, sep=" ")

#rename exonerate result columns
colnames(exonerate_rslt) <- c("query", "query_start", "query_end", "scaffold", "scaffold_start", 
                          "scaffold_end", "strand", "exonerate_score","intron_number","best_query", "evalue")



#Load blast results , rename columns, and extract best-hit regions
blast_rslt <- read.table("T1R_vs_Genome.blastn", header=FALSE, sep="\t")

colnames(blast_rslt) <- c("query", "sseqid", "pident", "length", "mismatch", 
                          "gapopen", "qstart", "qend", "sstart", "send",
                          "evalue", "bitscore")


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

blast_rslt <- blast_rslt %>% select(sseqid, true_start, true_end, evalue, strand, length)
colnames(blast_rslt) <- c("seqnames", "start", "end", "evalue", "strand", "length")

#keep only blast with an evalue < 1e-10 
blast_rslt <- blast_rslt %>% filter(evalue <= 1e-10) %>% mutate(scaffold_length = end-start) %>% filter(scaffold_length > 50)



blast_rslt_irange <- blast_rslt %>% as_granges()
blast_rslt_disjoin <- reduce(blast_rslt_irange,with.revmap=TRUE)
blast_result_besthits <- as.data.frame(blast_rslt_disjoin)

#blast_result_besthits <- blast_result_besthits %>% filter(width > 50) #keep only besthits longer than 80 nucleotides. Remove spurious blast results


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



exonerate_rslt <- exonerate_rslt %>% mutate(renamed_scaff = str_replace(scaffold, ":", "-"))  %>% separate(renamed_scaff, c("true_scaffold_name", "true_scaff_start", "true_scaff_stop"), "-") %>% mutate(true_coord_start = as.numeric(true_scaff_start)+scaffold_start) %>% mutate(true_coord_end = as.numeric(true_scaff_start)+scaffold_end)





#Remove exonerate results that have more exons than the number of blast hit detected in their region

colnames(exonerate_rslt) <- c("query", "query_start", "query_end", "false_scaffold_name", "false_scaffold_start", 
                          "false_scaffold_end", "strand", "exonerate_score","intron_number", "best_query", "evalue","seqnames", "true_scaff_start", "true_scaff_stop","start", "end")



#Remove genes with more blast hits that the number of exon, most probably chimeric genes

exonerate_rslt_filtered <- left_join(exonerate_rslt, blast_result_besthits, by = "seqnames", suffix=c("", "_b")) %>%
	filter(start_b >= start-50) %>%
	filter(end_b <= end+50)%>%
	group_by(seqnames, start, end, intron_number, query, query_end, query_start, evalue, best_query, strand, exonerate_score) %>%
	count(name="n_hit")%>%
	ungroup() %>% 
	filter(intron_number+1 >= n_hit) 


exonerate_rslt_filtered <- as.data.frame(exonerate_rslt_filtered %>% mutate(query_length = query_end - query_start) %>% mutate(scaffold_length = end-start) %>% filter(intron_number >= 1))




#Make filters on gene length and on scaffold length to find best hits
#exonerate_rslt_filtered %>% filter(seqnames=="groupI") %>% filter(start > 10738000) %>% filter(end < 10747784)

Best_hits_filtered <- as.data.frame(NULL) 

for(putative_protein_length in c(750,700,650,500, 400,200,100)){ 

	exonerate_rslt_filtered_round <- exonerate_rslt_filtered %>% filter(query_length > putative_protein_length)

	
	if(nrow(Best_hits_filtered) > 0){
		for (row in 1:nrow(Best_hits_filtered)) {
			
			
			exonerate_rslt_filtered_round <- exonerate_rslt_filtered_round %>% filter(!(seqnames == as.character(Best_hits_filtered[row, "seqnames"]) & start <= (	Best_hits_filtered[row, "start"]) & end >= (Best_hits_filtered[row, "start"])))
			
			
			exonerate_rslt_filtered_round <- exonerate_rslt_filtered_round %>% filter(!(seqnames == as.character(Best_hits_filtered[row, "seqnames"]) & start >= (	Best_hits_filtered[row, "start"]) & end <=(Best_hits_filtered[row, "end"])))
			
			exonerate_rslt_filtered_round <- exonerate_rslt_filtered_round %>% filter(!(seqnames == as.character(Best_hits_filtered[row, "seqnames"]) & start <= (	Best_hits_filtered[row, "end"]) & end >= (Best_hits_filtered[row, "end"])))
			
			exonerate_rslt_filtered_round <- exonerate_rslt_filtered_round %>% filter(!(seqnames == as.character(Best_hits_filtered[row, "seqnames"]) & start <= (	Best_hits_filtered[row, "start"]) & end >= (Best_hits_filtered[row, "end"])))
			
			  
		}
	}
	
		
	while (nrow(exonerate_rslt_filtered_round) > 0){
	
	
		exonerate_rslt_filtered_round_irange <- exonerate_rslt_filtered_round %>% as_granges()
		exonerate_rslt_filtered_round_disjoin <- reduce(exonerate_rslt_filtered_round_irange,with.revmap=TRUE)
		list_revmap <- as.data.frame(mcols(exonerate_rslt_filtered_round_disjoin))
	
		filtered_data <- c()
		for(i in 1:nrow(list_revmap)){
	  	filtered_data <- c(filtered_data, (slice(list_revmap, i) %>% unlist(use.names=FALSE))[which.min(slice(exonerate_rslt_filtered_round, slice(list_revmap, i) %>% unlist(use.names=FALSE))$scaffold_length)])
		}
	
	
		Best_hits_filtered <- rbind(Best_hits_filtered, slice(exonerate_rslt_filtered_round, filtered_data))
	
	
		if(nrow(Best_hits_filtered) > 0){
			for (row in 1:nrow(Best_hits_filtered)) {
			
			
				exonerate_rslt_filtered_round <- exonerate_rslt_filtered_round %>% filter(!(seqnames == as.character(Best_hits_filtered[row, "seqnames"]) & start <= (	Best_hits_filtered[row, "start"]) & end >= (Best_hits_filtered[row, "start"])))
			
			
				exonerate_rslt_filtered_round <- exonerate_rslt_filtered_round %>% filter(!(seqnames == as.character(Best_hits_filtered[row, "seqnames"]) & start >= (	Best_hits_filtered[row, "start"]) & end <=(Best_hits_filtered[row, "end"])))
			
				exonerate_rslt_filtered_round <- exonerate_rslt_filtered_round %>% filter(!(seqnames == as.character(Best_hits_filtered[row, "seqnames"]) & start <= (	Best_hits_filtered[row, "end"]) & end >= (Best_hits_filtered[row, "end"])))
			
				exonerate_rslt_filtered_round <- exonerate_rslt_filtered_round %>% filter(!(seqnames == as.character(Best_hits_filtered[row, "seqnames"]) & start <= (	Best_hits_filtered[row, "start"]) & end >= (Best_hits_filtered[row, "end"])))
			
			  
			}
		
		}
	
	
	}
	

}

#Best_hits_filtered_F %>% filter(seqnames == "NC_053154.1") %>% filter(start > 123737100) %>% filter(end < 123753045)
#exonerate_rslt_filtered %>% filter(seqnames == "NC_053154.1") %>% filter(start > 123737100) %>% filter(end < 123753045)


## Find the best query per exonerate best hits

if (nrow(Best_hits_filtered) > 0) {

	Best_hits_filtered <- Best_hits_filtered %>% mutate(extanded_start = start-1500) %>% mutate(extanded_end = end+1500)


	Best_hits_filtered <- Best_hits_filtered %>% mutate(new_coord_start =  case_when(
	  extanded_start <= 1 ~ 1,
	  extanded_start > 1 ~ extanded_start
	))

		
		
	Best_queries_list <- c()
	for (row in 1:nrow(Best_hits_filtered)) {
	
		current_scaffold <- Best_hits_filtered[row,]$seqnames
		start_region <- Best_hits_filtered[row,]$new_coord_start
		end_region <- Best_hits_filtered[row,]$extanded_end
		curr_query_length <- Best_hits_filtered[row,]$query_length
	
		best_query <- head(exonerate_rslt_filtered %>% filter(query_length >= curr_query_length) %>% filter(seqnames == current_scaffold) %>% filter(start >= start_region) %>% filter(end <= end_region) %>% arrange(evalue, desc(query_length)), 1)  %>% pull(best_query)
	
		Best_queries_list <- c(Best_queries_list, best_query)
	
	}

	Best_hits_filtered <- cbind(Best_hits_filtered, Best_queries_list)

	Best_hits_filtered_F <- Best_hits_filtered %>% dplyr::select(seqnames, new_coord_start, extanded_end, strand, exonerate_score, query, Best_queries_list)
	

	} else { 
	
	Best_hits_filtered_F <- as.data.frame(NULL) 
	
}




#Problem, check : AFNZ01064187.1


#write the result in a table
write.table(Best_hits_filtered_F, file="Parsed_exonerate_gene_regions.tsv", quote=FALSE, sep='\t', row.names = FALSE, col.names=FALSE)

#groupI:10738727-10747792




