#load packages
library(data.table)
library(dplyr)
library("plyranges")
library("GenomicRanges")




#load the file containing regions with genes to be extended
coordinates_to_extend <- read.table("Coordinates_genes_verifications.tsv", header=FALSE, sep="\t")
#rename columns
colnames(coordinates_to_extend) <- c("seqnames", "start", "end")


#load the file containing regions with correct genes
regions_correct_genes <- read.table("Coordinates_correct_genes.tsv", header=FALSE, sep="\t")
#rename columns
colnames(regions_correct_genes) <- c("seqnames", "start", "end")




new_coord_start_v <- c()
new_coord_end_v <- c()

for (row in 1:nrow(coordinates_to_extend)) {
  scaffold_name <- coordinates_to_extend[row, "seqnames"]
  curr_start  <-  coordinates_to_extend[row, "start"]
  curr_end <-  coordinates_to_extend[row, "end"]
  
  if (nrow(regions_correct_genes %>% filter(seqnames == scaffold_name) %>% filter(end >= (curr_start-25000) & start < curr_end)) > 0){

    new_curr_start = (tail(regions_correct_genes %>% filter(seqnames == scaffold_name) %>% filter(end >= (curr_start-25000) & start < curr_start) %>% dplyr::arrange(end), 1)$end) + 100

  } else new_curr_start = curr_start-25000



  if (nrow(regions_correct_genes %>% filter(seqnames == scaffold_name) %>% filter(start <= (curr_end+25000) & end > curr_end)) > 0){

    new_curr_end = (head(regions_correct_genes %>% filter(seqnames == scaffold_name) %>% filter(start <= (curr_end+25000) & end > curr_end) %>% dplyr::arrange(start), 1)$start) - 100

  } else new_curr_end = curr_end+25000



  new_coord_start_v <- c(new_coord_start_v, new_curr_start)
  new_coord_end_v <- c(new_coord_end_v, new_curr_end)

}




coordinates_to_extend <- coordinates_to_extend %>% mutate(Extanded_start = new_coord_start_v)
coordinates_to_extend <- coordinates_to_extend %>% mutate(Extanded_end = new_coord_end_v)


coordinates_to_extend_test <- coordinates_to_extend %>% select(seqnames, Extanded_start, Extanded_end)
colnames(coordinates_to_extend_test) <- c("seqnames", "start", "end")


coordinates_to_extend_test_irange <- coordinates_to_extend_test %>% as_granges()

#reduce the table
coordinates_to_extend_test_disjoin <- reduce(coordinates_to_extend_test_irange,with.revmap=TRUE)



coordinates_to_extend_test_disjoin_df <- as.data.frame(coordinates_to_extend_test_disjoin)



coordinates_to_extend_test_disjoin_df <- coordinates_to_extend_test_disjoin_df %>% mutate(new_coord_start = case_when(
  as.numeric(start) <= 0 ~ 1,
  as.numeric(start) > 0 ~ as.numeric(start)))



coordinates_to_extend_test_disjoin_df$samtools_name <- paste(coordinates_to_extend_test_disjoin_df$seqnames, ":", coordinates_to_extend_test_disjoin_df$new_coord_start, "-", coordinates_to_extend_test_disjoin_df$end, sep = "")


#Write regions in a text file
write(coordinates_to_extend_test_disjoin_df$samtools_name, file="Extended_regions.tsv")
