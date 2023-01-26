#load packages

library("ape")
library("dplyr")
library("phytools")

#load the tree
mytree <- read.tree("Putative_or_plus_known_or_plus_outgroup.prot.aln.treefile")

#name internal nodes of the tree
mytree_label <- makeNodeLabel(mytree, method = "number", prefix="Node")

#import the name of known OR genes
known_OR_genes_typeI <- scan("Known_OR_genes_TYPEI.id", what="character")
known_OR_genes_typeII <- scan("Known_OR_genes_TYPEII.id", what="character")
known_OR_genes_all <- c(known_OR_genes_typeI, known_OR_genes_typeII)


#Root the tree at the common ancestor of every known OR genes
MRCA_OR <- findMRCA(mytree_label, tips=known_OR_genes_all, type="node") #extract the node name
mytree_rooted <- root(mytree_label, node=MRCA_OR, resolve.root= TRUE) #root
MRCA_OR <- findMRCA(mytree_rooted, tips=known_OR_genes_all, type="node") #extract the name of the ancestral node to all OR genes

#extract all OR genes
OR_genes <- as.vector(extract.clade(mytree_rooted, MRCA_OR)$tip)

#Remove from the OR gene list the name of already known OR genes
Current_species_OR <- setdiff(OR_genes, known_OR_genes_all)

#write in file
write(x=Current_species_OR, file="Current_species_OR.txt")



