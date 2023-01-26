library("ape")
library("dplyr")
library("phytools")

#Test pipeline tree refining

mytree <- read.tree("Putative_or_plus_known_or_plus_outgroup.prot.aln.treefile")

mytree_label <- makeNodeLabel(mytree, method = "number", prefix="Node")

known_OR_genes_typeI <- scan("Known_OR_genes_TYPEI.id", what="character")
known_OR_genes_typeII <- scan("Known_OR_genes_TYPEII.id", what="character")

known_OR_genes_all <- c(known_OR_genes_typeI, known_OR_genes_typeII)

#Root the tree at the common ancestor of every OR genes

MRCA_ORs_typeI <- findMRCA(mytree_label, tips=known_OR_genes_typeI, type="node")
mytree_rooted <- root(mytree_label, node=MRCA_ORs_typeI, resolve.root= TRUE)
MRCA_ORs_typeI <- findMRCA(mytree_rooted, tips=known_OR_genes_typeI, type="node")
tpyeI_genes <- extract.clade(mytree_rooted, node=MRCA_ORs_typeI)


MRCA_ORs_typeII <- findMRCA(mytree_label, tips=known_OR_genes_typeII, type="node")
mytree_rooted <- root(mytree_label, node=MRCA_ORs_typeII, resolve.root= TRUE)
MRCA_ORs_typeII <- findMRCA(mytree_rooted, tips=known_OR_genes_typeII, type="node")
typeII_genes <- extract.clade(mytree_rooted, node=MRCA_ORs_typeII)


OR_genes <- c(tpyeI_genes$tip.label, typeII_genes$tip.label)

Current_species_OR <- setdiff(OR_genes, known_OR_genes_all)

write(x=Current_species_OR, file="Current_species_OR.txt")


