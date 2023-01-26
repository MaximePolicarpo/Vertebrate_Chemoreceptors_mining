#load packages
library(data.table)
library(dplyr)
library("plyranges")
library("GenomicRanges")


#load coordinates of found genes 

functionnal_genes <- read.table("Coordinates_genes_already_found.tsv", header=FALSE, sep="\t")
colnames(functionnal_genes) <- c("seqnames", "true_start", "true_end")

#re-orientate results and put an identifier to genes
functionnal_genes <- functionnal_genes %>% mutate(start = case_when(
  true_start < true_end ~ true_start,
  true_end < true_start ~ true_end
))
functionnal_genes <- functionnal_genes %>% mutate(end = case_when(
  true_start < true_end ~ true_end,
  true_end < true_start ~ true_start
))
functionnal_genes <- functionnal_genes %>% mutate(gene_state = "F")



#Merge the data tables

functionnal_genes <- functionnal_genes %>% mutate(length = end - start)


#Put the table as a grange object
all_genes_df_irange <- functionnal_genes %>% as_granges()


#Reduce the table to merge overlapping results
all_genes_df_disjoin <- reduce(all_genes_df_irange,with.revmap=TRUE)



#for each overlapping regions, keep only the longest gene


list_revmap <- as.data.frame(mcols(all_genes_df_disjoin))



filtered_data <- c()
for(i in 1:nrow(list_revmap)){
  filtered_data <- c(filtered_data, (slice(list_revmap, i) %>% unlist(use.names=FALSE))[which.max(slice(functionnal_genes,  slice(list_revmap, i) %>% unlist(use.names=FALSE))$length)])
}


best_genes_df <- slice(functionnal_genes, filtered_data)



#Seperate table depending on gene state and export files

best_genes_functionnal <- best_genes_df 

write.table(best_genes_functionnal, file="Best_combination_functional.tsv", quote=FALSE, sep='\t', row.names = FALSE, col.names=FALSE)


