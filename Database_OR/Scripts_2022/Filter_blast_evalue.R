library(dplyr)

blast_rslt <- read.table("Initial_OR_vs_Genome.blastn", header=FALSE, sep="\t")

#rename blast result columns
colnames(blast_rslt) <- c("query", "sseqid", "pident", "length", "mismatch", 
                          "gapopen", "qstart", "qend", "sstart", "send",
                          "evalue", "bitscore", "qframe", "sframe", "qlen", "slen")


#filter blast results
blast_rslt_filtered <- blast_rslt %>% filter(evalue <= 1e-20)

#filter columns
blast_rslt_normal <- blast_rslt %>% dplyr::select(query, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore)

#write tables
write.table(blast_rslt_filtered, file="tblastn_functionnal_or_vs_genome.blastn", quote=FALSE, sep='\t', row.names = FALSE, col.names=FALSE)
write.table(blast_rslt_normal, file="OR_vs_Genome.blastn", quote=FALSE, sep='\t', row.names = FALSE, col.names=FALSE)

