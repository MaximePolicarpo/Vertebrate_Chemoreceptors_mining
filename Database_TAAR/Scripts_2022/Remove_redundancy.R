#load packages
library(data.table)
library(dplyr)
library("plyranges")
library("GenomicRanges")


#load coordinates of found genes 

if (file.exists("Coordinates_genes_final.tsv")){
	functionnal_genes <- read.table("Coordinates_genes_final.tsv", header=FALSE, sep="\t")
} else functionnal_genes <- as.data.frame(matrix(ncol = 3, nrow = 0))


if (file.exists("Coordinates_ambigous_final.tsv")){
	ambigous_genes <- read.table("Coordinates_ambigous_final.tsv", header=FALSE, sep="\t")
} else ambigous_genes <- as.data.frame(matrix(ncol = 3, nrow = 0))


if (file.exists("Coordinates_pseudogenes_final.tsv")){
	pseudo_genes <- read.table("Coordinates_pseudogenes_final.tsv", header=FALSE, sep="\t")
} else pseudo_genes <- as.data.frame(matrix(ncol = 3, nrow = 0))

#rename columns

colnames(functionnal_genes) <- c("seqnames", "true_start", "true_end")
colnames(ambigous_genes) <- c("seqnames", "true_start", "true_end")
colnames(pseudo_genes) <- c("seqnames", "true_start", "true_end")





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





ambigous_genes <- ambigous_genes %>% mutate(start = case_when(
  true_start < true_end ~ true_start,
  true_end < true_start ~ true_end
))
ambigous_genes <- ambigous_genes %>% mutate(end = case_when(
  true_start < true_end ~ true_end,
  true_end < true_start ~ true_start
))
ambigous_genes <- ambigous_genes %>% mutate(gene_state = "A")




pseudo_genes <- pseudo_genes %>% mutate(start = case_when(
  true_start < true_end ~ true_start,
  true_end < true_start ~ true_end
))
pseudo_genes <- pseudo_genes %>% mutate(end = case_when(
  true_start < true_end ~ true_end,
  true_end < true_start ~ true_start
))
pseudo_genes <- pseudo_genes %>% mutate(gene_state = "P")




#Merge the data tables

all_genes_df <- do.call("rbind", list(functionnal_genes, pseudo_genes, ambigous_genes))


all_genes_df <- all_genes_df %>% mutate(length = end - start)


#Put the table as a grange object
all_genes_df_irange <- all_genes_df %>% as_granges()


#Reduce the table to merge overlapping results
all_genes_df_disjoin <- reduce(all_genes_df_irange,with.revmap=TRUE)



#for each overlapping regions, keep only the longest gene


list_revmap <- as.data.frame(mcols(all_genes_df_disjoin))



filtered_data <- c()
for(i in 1:nrow(list_revmap)){
  filtered_data <- c(filtered_data, (slice(list_revmap, i) %>% unlist(use.names=FALSE))[which.max(slice(all_genes_df,  slice(list_revmap, i) %>% unlist(use.names=FALSE))$length)])
}


best_genes_df <- slice(all_genes_df, filtered_data)



#Seperate table depending on gene state and export files

best_genes_functionnal <- best_genes_df %>% filter(gene_state == "F")
best_genes_ambigous <- best_genes_df %>% filter(gene_state == "A")
best_genes_pseudogenes <- best_genes_df %>% filter(gene_state == "P")



write.table(best_genes_functionnal, file="best_genes_functionnal.tsv", quote=FALSE, sep='\t', row.names = FALSE, col.names=FALSE)
write.table(best_genes_ambigous, file="best_genes_ambigous.tsv", quote=FALSE, sep='\t', row.names = FALSE, col.names=FALSE)
write.table(best_genes_pseudogenes, file="best_genes_pseudogenes.tsv", quote=FALSE, sep='\t', row.names = FALSE, col.names=FALSE)




