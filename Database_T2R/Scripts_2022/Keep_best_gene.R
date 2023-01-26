#load packages
library(data.table)
library(dplyr)
library("plyranges")
library("GenomicRanges")




#load the file containing genes and their coordinates
gene_coordinates <- read.table("Combined_old_new_genes_to_parse.tsv", header=FALSE, sep="\t")
#rename columns
colnames(gene_coordinates) <- c("seqnames", "start", "end", "exon_number","gene_name")



gene_coordinates_irange <- gene_coordinates %>% as_granges()

#reduce the table
gene_coordinates_disjoin <- reduce(gene_coordinates_irange,with.revmap=TRUE)

list_revmap <- as.data.frame(mcols(gene_coordinates_disjoin))


#keep genes with the maximum exon number

filtered_data <- c()
for(i in 1:nrow(list_revmap)){
  filtered_data <- c(filtered_data, (slice(list_revmap, i) %>% unlist(use.names=FALSE))[which.max(slice(gene_coordinates,  slice(list_revmap, i) %>% unlist(use.names=FALSE))$exon_number)])
}

Best_genes <- slice(gene_coordinates, filtered_data)



write(Best_genes$gene_name, file="Genes_to_keep.tsv")
