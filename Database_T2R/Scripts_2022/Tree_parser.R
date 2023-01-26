#load packages

library("ape")
library("dplyr")
library("phytools")

#load the tree
mytree <- read.tree("Putative_t2r_plus_known_t2r_plus_outgroup.prot.aln.treefile")

#name internal nodes of the tree
mytree_label <- makeNodeLabel(mytree, method = "number", prefix="Node")

#import the name of known T2R genes
known_T2R_genes <- scan("T2R_genes_tree.id", what="character")
known_V1R_genes <- scan("V1R_genes_tree.id", what="character")

#Root the tree at the common ancestor of outgroup sequences
MRCA_outgroup <- findMRCA(mytree_label, tips=known_V1R_genes, type="node")
mytree_rooted <- root(mytree_label, node=MRCA_outgroup, resolve.root= TRUE) #root

#Check the MRCA or V2R genes

MRCA_T2R <- findMRCA(mytree_rooted, tips=known_T2R_genes, type="node")

#grep all tips from these MRCS
T2R_genes <- as.vector(extract.clade(mytree_rooted, MRCA_T2R)$tip)

#Remove alrdy known TAAR genes
Current_species_T2R <- setdiff(T2R_genes, known_T2R_genes)


write(x=Current_species_T2R, file="Current_species_T2R.txt")

