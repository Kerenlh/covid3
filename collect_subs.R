setwd("/Users/keren/Dropbox/covid")
source("COVID_functions.R")

setwd("/Users/keren/Dropbox/covid/vars")
aligned_seqs = readDNAStringSet(file="aligned_tree_seqs_decipher")
tree = read.tree("tree.nwk")
tree = short_names_tree(tree)

# contrast = matrix(data = c(1,0,0,0,0,
#                            0,1,0,0,0,
#                            0,0,1,0,0,
#                            0,0,0,1,0,
#                            1,0,1,0,0,
#                            0,1,0,1,0,
#                            0,0,0,0,1,
#                            1,1,1,1,0,
#                            1,1,1,1,1),
#                   ncol = 5, byrow = TRUE)
# dimnames(contrast) = list(c("A","C","G","T","R","Y","-","N","?"),
#                           c("A", "C", "G", "T", "-"))

aligned_seqs = as.DNAbin(aligned_seqs)
phyDat_seqs = as.phyDat(aligned_seqs)
# phyDat_seqs = as.phyDat(aligned_seqs,type = "USER",contrast = contrast,levels=c("A","C","G","T","-"),
#                         ambiguity = c("?", "N"))
phyDat_seqs = short_names_phyDat(phyDat_seqs)


site_patterns = find_site_patterns(aligned_seqs)
# save(site_patterns, file = "site_patterns")
load("site_patterns")

all_site_patterns = unique(site_patterns)
site_pattern_indices = matrix(0,length(site_patterns),1)
for (i in 1:length(all_site_patterns)){
  site_pattern_indices[which(site_patterns==all_site_patterns[i])] = i
}

phyDat_patterns = NULL
for (i in 1:length(phyDat_seqs[[1]])){
  phyDat_pattern = NULL
  for (j in 1:length(phyDat_seqs)){
    phyDat_pattern = c(phyDat_pattern,phyDat_seqs[[j]][i])
  }
  phyDat_patterns = rbind(phyDat_patterns,phyDat_pattern)
}

reconstructed_seqs = fitch(tree = tree, data = phyDat_seqs,site = "data")

#collect_subs
subs = NULL
sites = NULL
for (i in 1:dim(tree$edge)[1]){
  # if (tree$edge[i,2]<=length(tree$tip.label)){
  #   rel_accesion = tree$tip.label[tree$edge[i,2]]
  # }else{
  #   rel_accesion = NA
  # }
  subs_plcs = reconstructed_seqs$dat[,tree$edge[i,1]]!=reconstructed_seqs$dat[,tree$edge[i,2]]
  subs_plcs_plcs = which(subs_plcs==TRUE)
  curr_sites = NULL
  num_sites = NULL
  for (j in 1:length(subs_plcs_plcs)){
    curr_sites = rbind(curr_sites,list(which(site_patern_indices==subs_plcs_plcs[j])))
    num_sites = rbind(num_sites, length(which(site_patern_indices==subs_plcs_plcs[j])))
  }
  sites = rbind(sites,curr_sites)
  subs = rbind(subs, cbind(which(subs_plcs==TRUE),reconstructed_seqs$dat[subs_plcs,tree$edge[i,1]],
                           reconstructed_seqs$dat[subs_plcs,tree$edge[i,2]], 
                           matrix(tree$edge[i,1],sum(subs_plcs),1),
                           matrix(tree$edge[i,2],sum(subs_plcs),1) ,num_sites) )
}
colnames(subs) = c("site_pattern","before","after","f_edge","s_edge","num_sites")
subs = data.frame(subs)

subs$before = change_to_letters(subs$before)
subs$after = change_to_letters(subs$after)

subs_muts = apply(cbind(subs$before,matrix("_",dim(subs)[1],1),subs$after), 1, paste0,
                  collapse="")  