#load packages

library("ape")
library("dplyr")
library("phytools")

#load the tree
mytree <- read.tree("Pseudogene_verif_alignment.aln.treefile")

#name internal nodes of the tree
mytree_label <- makeNodeLabel(mytree, method = "number", prefix="Node")

#import the name of known OR genes
known_TAAR_L <- scan("taar_l.id", what="character")
known_TAAR <- scan("TAAR_known.id", what="character")
Outgroups_seqs <- scan("Outgroups_id.txt", what="character")


#Root the tree at the common ancestor of outgroup sequences
MRCA_outgroup <- findMRCA(mytree_label, tips=Outgroups_seqs, type="node")
mytree_rooted <- root(mytree_label, node=MRCA_outgroup, resolve.root= TRUE) #root

#Check the MRCA or TAAR-genes and TAAR_like genes

MRCA_TAAR <- findMRCA(mytree_rooted, tips=known_TAAR, type="node")
MRCA_TAAR_L <- findMRCA(mytree_rooted, tips=known_TAAR_L, type="node")

#grep all tips from these MRCS
TAAR_genes <- as.vector(extract.clade(mytree_rooted, MRCA_TAAR)$tip)
TAAR_L_genes <- as.vector(extract.clade(mytree_rooted, MRCA_TAAR_L)$tip)


#combine data
All_TAAR <- c(TAAR_genes,TAAR_L_genes)

#Remove alrdy known TAAR genes
known_taar_all <- c(known_TAAR, known_TAAR_L)
Current_species_TAAR <- setdiff(All_TAAR, known_taar_all)


write(x=Current_species_TAAR, file="Current_species_pseudo_TAAR.txt")





#MRCA_OR <- findMRCA(mytree_label, tips=known_OR_genes_all, type="node") #extract the node name
#mytree_rooted <- root(mytree_label, node=MRCA_OR, resolve.root= TRUE) #root
#MRCA_OR <- findMRCA(mytree_rooted, tips=known_OR_genes_all, type="node") #extract the name of the ancestral node to all OR genes

#extract all OR genes
#OR_genes <- as.vector(extract.clade(mytree_rooted, MRCA_OR)$tip)

#Remove from the OR gene list the name of already known OR genes
#Current_species_OR <- setdiff(OR_genes, known_OR_genes_all)

#write in file
#write(x=Current_species_OR, file="Current_species_OR.txt")



