setwd("/Users/keren/Dropbox/EC2_files")
source("COVID_functions.R")
aligned_seqs = readDNAStringSet(file="aligned_tree_seqs_decipher")
setwd("/Users/keren/Dropbox/covid")
tree = read.tree("aligned_tree_seqs_decipher tree.newick")
tree = short_names_tree(tree)

# setwd("/Users/keren/Dropbox/covid/vars/missing_site_patterns2/")
# load("missing_all_site_patterns")
setClass(Class='Tree', 
         representation(name= 'character',edge = 'numeric', parent_edge = 'numeric',
                        parent_name = 'character',  
                        children_edges = 'matrix', children_names = 'character', 
                        probs = 'matrix', updated_probs = 'matrix', base = 'character',leaf = 'logical'))


find_details_not_leaf = function(curr_name,curr_site_pattern,
                                 tree,tree_list){
  #print(curr_name)
  tmp = strsplit(curr_name,split = "")
  if (paste0(tmp[[1]][1:8],collapse = "")=="no_name_"){
    curr_edge = as.numeric(paste0(tmp[[1]][9:length(tmp[[1]])],collapse = ""))
  }else{
    curr_edge = which(tree$tip.label==curr_name)
    leaf_bug_check = 1
  }
  curr_children_edges = as.matrix(tree$edge[which(tree$edge[,1]==curr_edge),2])
  if (length(curr_children_edges)==0){ # This is a leaf
    if (leaf_bug_check==0){
      print("This is a leaf identification bug")
    }
    tree_list = find_details_leaf(curr_name,curr_site_pattern,tree,tree_list)
    return(tree_list)
  }else{
    leaf = FALSE
    curr_parent_edge = tree$edge[which(tree$edge[,2]==curr_edge),1]
    curr_parent_name = tree$tip.label[curr_parent_edge]
    if (curr_edge == 12225){ # This is the root
      curr_parent_name = "no_name_0"
    }
    if(is.na(curr_parent_name)){
      curr_parent_name = paste0("no_name_",curr_parent_edge,collapse = "")
    }
    curr_children_names = NULL
    for (i in 1:length(curr_children_edges)){
      if (curr_children_edges[i]<=length(tree$tip.label)){
        curr_children_names = c(curr_children_names,tree$tip.label[curr_children_edges[i]])
      }else{
        curr_children_names = c(curr_children_names,paste0("no_name_",curr_children_edges[i]))
      }
    }
    all_probs = matrix(0,1,5)
    for (i in 1:length(curr_children_names)){
      if (is.null(tree_list[[curr_children_names[i]]])){
        tree_list = find_details_not_leaf(curr_children_names[i],curr_site_pattern,
                                          tree,tree_list)
      } 
      all_probs = tree_list[[curr_children_names[i]]]@probs + all_probs
    }
    probs = all_probs/length(curr_children_names)
    updated_probs = matrix(0,1,5)
    base = colnames(probs)[which.max(probs)]
    if(length(which(probs==max(probs)))>1){
      print(paste(curr_name, probs))
    }
    tt = new("Tree", name = curr_name, edge = curr_edge, parent_edge = curr_parent_edge,
             parent_name = curr_parent_name, children_edges = curr_children_edges,
             children_names = curr_children_names, base = base, probs = probs,
             updated_probs = updated_probs, leaf = leaf)
    tree_list[[curr_name]] = tt
    return(tree_list)
  }
}

find_details_leaf = function(curr_name,curr_site_pattern,tree,
                             tree_list){
  print(curr_name)
  curr_edge = which(tree$tip.label==curr_name)
  curr_parent_edge = tree$edge[which(tree$edge[,2]==curr_edge),1]
  curr_parent_name = tree$tip.label[curr_parent_edge]
  if(is.na(curr_parent_name)){
    curr_parent_name = paste0("no_name_",curr_parent_edge,collapse = "")
  }
  curr_children_edges = tree$edge[which(tree$edge[,1]==curr_edge),2]
  if (length(curr_children_edges)>0 | curr_edge>length(tree$tip.label)){
    print("bug, this is not a leaf!")
  }else{
    base = "A" 
    # base = curr_site_pattern[which(seqs_names==curr_name)]
    probs = get_probs_from_base(base)
    updated_probs = matrix(0,1,5)
    leaf = TRUE
  }
  tt = new("Tree", name = curr_name, edge = curr_edge, parent_edge = curr_parent_edge,
           parent_name = curr_parent_name, children_edges = matrix(0,1,1),
           children_names = "_", base = base, probs = probs,
           updated_probs = updated_probs, leaf = leaf)
  tree_list[[curr_name]] = tt
  return(tree_list)
}

get_probs_from_base = function(base){
  probs = matrix(0,1,5)
  colnames(probs) = c("A","C","G","T","-")
  if (is.element(base,colnames(probs))){
    probs[,base] = 1
  }else if (base=="R"){
    probs[,c("A","G")] = 0.5
  }else if (base=="Y"){
    probs[,c("C","T")] = 0.5
  }else if (base=="S"){
    probs[,c("C","G")] = 0.5
  }else if (base=="W"){
    probs[,c("A","T")] = 0.5
  }else if (base=="K"){
    probs[,c("T","G")] = 0.5
  }else if (base=="M"){
    probs[,c("A","C")] = 0.5
  }else if (base=="B"){
    probs[,c("C","T","G")] = 1/3
  }else if (base=="D"){
    probs[,c("A","T","G")] = 1/3
  }else if (base=="H"){
    probs[,c("C","T","A")] = 1/3
  }else if (base=="V"){
    probs[,c("C","A","G")] = 1/3
  }else if (base=="N"){
    probs[,c("C","A","G","T")] = 1/4
  }
  return(probs)
}

update_probs = function(tree_list,curr_name,epsilon){
  if (tree_list[[curr_name]]@leaf==FALSE){
    children = tree_list[[curr_name]]@children_names
    for (i in 1:length(children)){
      tree_list[[children[i]]]@updated_probs = 
        (tree_list[[children[i]]]@probs + epsilon*tree_list[[curr_name]]@probs)/
        (1+epsilon)
      tmp = tree_list[[children[i]]]@probs
      if(length(which(tmp==max(tmp)))>1){
        print(paste(children[i], tmp, tree_list[[children[i]]]@updated_probs))
        tree_list[[children[i]]]@base = 
          colnames(tmp)[which.max(tree_list[[children[i]]]@updated_probs)]
        print(tree_list[[children[i]]]@base)
      }
      tree_list = update_probs(tree_list,children[i],epsilon)        
    }
  }
  return(tree_list) 
}

seqs_names = names(aligned_seqs)
root_edge = unique(tree$edge[which(is.element(tree$edge[,1],tree$edge[,2])==FALSE),1])
#root_edge = 12225
root_name = paste0(c("no_name_",root_edge),collapse = "")
# setwd("/Users/keren/Dropbox/covid/vars/missing_site_patterns2/trees")
# setwd("/home/ubuntu/Dropbox/EC2_files/trees")
# unique_site_patterns = unique(all_site_patterns)
# site_patterns_numbers = c(13351:13366)
# plcs = 1:length(unique_site_patterns)
epsilon = 0.01
# for (j in plcs[((parallel_nums-1)*cluster_delta+1)]:plcs[parallel_nums*cluster_delta]){#  1:length(unique_site_patterns)){
for (j in 1:length(unique_site_patterns)){
  print(j)
  tree_list = list()
  curr_site_pattern = strsplit(unique_site_patterns[j],split = "")[[1]]
  start_time = Sys.time()
  tree_list = find_details_not_leaf(root_name,curr_site_pattern,
                                    tree,tree_list)
  tree_list = update_probs(tree_list, root_name,epsilon)
  save(tree_list,file = paste0("tree_list_",site_patterns_numbers[j],collapse=""))
  end_time = Sys.time()
  print(end_time-start_time)
}
