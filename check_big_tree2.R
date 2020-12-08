setwd("/Users/keren/Dropbox/covid/vars/big_trees")
load("all_big_trees3")

###########3
# Check diffreneces between seqs:
leaf_names = NULL
split_leafs = list()
for (i in 1:length(all_tree_list)){
  print(i)
  if (all_tree_list[[i]]@leaf==TRUE){
    leaf_names = rbind(leaf_names,c(i,all_tree_list[[i]]@name))
  }
}

leaf_plcs = as.numeric(leaf_names[,1])
split_leafs = list()
for (i in 1:length(leaf_plcs)){
  print(i)
  split_leafs[[i]] = strsplit(all_tree_list[[leaf_plcs[i]]]@probs,split = "")[[1]]
}

setwd("/Users/keren/Dropbox/covid/vars/")
save(split_leafs,file = "split_leafs")
save(leaf_plcs,file = "leaf_plcs")

diffs_count = matrix(-1,length(leaf_plcs),length(leaf_plcs))
for (i in 1:length(split_leafs)){
  print(i)
  start_time = Sys.time()
  for (j in (i+1):length(split_leafs)){
    diffs_count[i,j] = length(which(split_leafs[[i]]!=split_leafs[[j]]))
  }
  print(Sys.time()-start_time)
}


