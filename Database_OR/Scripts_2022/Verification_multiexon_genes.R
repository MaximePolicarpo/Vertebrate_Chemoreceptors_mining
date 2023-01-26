#load packages
library(data.table)
library(dplyr)
library("plyranges")
library("GenomicRanges")

#load and rename columns of blast result

blast_rslt <- as.data.frame(NULL)
blast_rslt <- tryCatch(read.table("tblastn_verif_multiexon.tsv", header = FALSE, sep = '\t'), error=function(e) NULL)


if (length(blast_rslt) > 0) {  


	#blast_rslt <- read.table("tblastn_verif_multiexon.tsv", header=FALSE, sep="\t")
	colnames(blast_rslt) <- c("query", "sseqid", "pident", "length", "mismatch", 
	                          "gapopen", "qstart", "qend", "sstart", "send",
	                          "evalue", "bitscore", "qframe", "sframe", "qlen", "slen")
	
	
	
	#re-orientate blast results and print the strand
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
	
	
	#rename columns
	
	
	colnames(blast_rslt) <- c("query", "seqnames", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "seq_start", "seq_end", "evalue", "bitscore", "qframe", "sframe", "qlen", "slen", "start", "end", "strand")
	
	
	#transform the table as irange object
	
	blast_rslt_irange <- blast_rslt %>% as_granges()
	
	#reduce the table
	blast_rslt_disjoin <- reduce(blast_rslt_irange,with.revmap=TRUE)
	
	
	blast_rslt_disjoin_df <- as.data.frame(blast_rslt_disjoin)
	
	blast_rslt_disjoin_df <- blast_rslt_disjoin_df %>% filter(width > 50)

	blast_rslt_disjoin_df <- blast_rslt_disjoin_df %>% select(seqnames, start, end)

} else {

	blast_rslt_disjoin_df <- as.data.frame(NULL)

}




write.table(blast_rslt_disjoin_df, file='blast_rslt_disjoin_df.tsv', quote=FALSE, sep='\t', row.names = FALSE,
            col.names=FALSE)


