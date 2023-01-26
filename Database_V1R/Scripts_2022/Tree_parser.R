#load packages

library("ape")
library("dplyr")
library("phytools")

#load the tree
mytree <- read.tree("Final_ALL_verification_alignment.aln.treefile")

#name internal nodes of the tree
mytree_label <- makeNodeLabel(mytree, method = "number", prefix="Node")

#import the name of known V1R genes

known_V1R <- scan("V1R_sequences.id", what="character")
Outgroups_seqs <- scan("T2R_id.txt", what="character")

#Root the tree at the common ancestor of outgroup sequences
MRCA_outgroup <- findMRCA(mytree_label, tips=Outgroups_seqs, type="node")
mytree_rooted <- root(mytree_label, node=MRCA_outgroup, resolve.root= TRUE) #root

#Check the MRCA or V1R genes

MRCA_V1R <- findMRCA(mytree_rooted, tips=known_V1R, type="node")

#grep all tips from these MRCS
V1R_genes <- as.vector(extract.clade(mytree_rooted, MRCA_V1R)$tip)

#Remove alrdy known TAAR genes
Current_species_V1R <- setdiff(V1R_genes, known_V1R)


write(x=Current_species_V1R, file="Current_species_V1R.txt")
