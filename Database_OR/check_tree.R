library("ape")
library("phytools")
args = commandArgs(trailingOnly=TRUE)
alpha_id <- args[1]
beta_id <- args[2]
delta_id <- args[3]
epsilon_id <- args[4]
eta_id <- args[5]
gamma_id <- args[6]
zeta_id <- args[7]
tree_file <- args[8]
output_file <- args[9]

mytree <- read.tree(tree_file)
outgroup_id <- scan("outgroup.id", "character")
outgroup_MRCA <- findMRCA(mytree, tips=outgroup_id, type="node")
mytree_rooted <- root(mytree, node=outgroup_MRCA, resolve.root= TRUE)

alpha <- scan(alpha_id , what="character")
beta <- scan(beta_id , what="character")
delta <- scan(delta_id , what="character")
epsilon <- scan(epsilon_id , what="character")
eta <- scan(eta_id , what="character")
gamma <- scan(gamma_id , what="character")
zeta <- scan(zeta_id , what="character")
all_to_prune <- c(alpha, beta, delta, epsilon, eta, gamma, zeta, outgroup_id)

mytree_rooted_pruned <- keep.tip(mytree_rooted, all_to_prune)

alpha_state=is.monophyletic(mytree_rooted_pruned, alpha)
beta_state=is.monophyletic(mytree_rooted_pruned, beta)
delta_state=is.monophyletic(mytree_rooted_pruned, delta)
epsilon_state=is.monophyletic(mytree_rooted_pruned, epsilon)
eta_state=is.monophyletic(mytree_rooted_pruned, eta)
gamma_state=is.monophyletic(mytree_rooted_pruned, gamma)
zeta_state=is.monophyletic(mytree_rooted_pruned, zeta)

all_states <- c(alpha_state, beta_state, delta_state, epsilon_state, eta_state, gamma_state, zeta_state)

print(all_states)
