setwd("/Users/keren/Dropbox/covid")
# setwd("C:/Users/ifog/Dropbox/covid")
source("COVID_functions.R")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DECIPHER")
library(DECIPHER)

tree = read.tree("tree.nwk")
tree = short_names_tree(tree)

seqs = readDNAStringSet(file = "sequences.fasta")
seqs = short_names_phyDat(seqs)

rel_seqs = return_only_tree_seqs(seqs,tree)
# aligned_seqs = msa(rel_seqs)

aligned <- AlignSeqs(rel_seqs)
setwd("/Users/keren/Dropbox/covid/vars")
# save(aligned,file = "aligned_tree_seqs_decipher")
# # view the alignment in a browser (optional)
BrowseSeqs(aligned, highlight=0)
# 
# # write the alignment to a new FASTA file
writeXStringSet(aligned,
                file="aligned_tree_seqs_decipher")


