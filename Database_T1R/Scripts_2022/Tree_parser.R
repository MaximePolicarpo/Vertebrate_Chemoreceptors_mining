#load packages

library("ape")
library("dplyr")
library("phytools")

#load the tree
mytree <- read.tree("Final_ALL_verification_alignment.aln.treefile")

#name internal nodes of the tree
mytree_label <- makeNodeLabel(mytree, method = "number", prefix="Node")

#import the name of known OR genes
known_T1R <- scan("Known_T1R_id.txt", what="character")
Outgroups_seqs <- scan("Known_outgroups_id.txt", what="character")


#Root the tree at the common ancestor of outgroup sequences
MRCA_outgroup <- findMRCA(mytree_label, tips=Outgroups_seqs, type="node")
mytree_rooted <- root(mytree_label, node=MRCA_outgroup, resolve.root= TRUE) #root

#Check the MRCA or T1R genes

MRCA_T1R <- findMRCA(mytree_rooted, tips=known_T1R, type="node")

#grep all tips from these MRCS
T1R_genes <- as.vector(extract.clade(mytree_rooted, MRCA_T1R)$tip)

#Remove alrdy known TAAR genes
Current_species_T1R <- setdiff(T1R_genes, known_T1R)


write(x=Current_species_T1R, file="Current_species_T1R.txt")








