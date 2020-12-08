setwd("/Users/keren/Dropbox/EC2_files/")
aligned_seqs = readDNAStringSet(file="aligned_tree_seqs_decipher")
tree = read.tree("tree.nwk")
setwd("/Users/keren/Dropbox/covid/vars/")
load("site_details")
load("missing_all_no_subs_plcs")
setwd("/Users/keren/Desktop/covid_files/missing_trees/")
missing_tree_names = dir()
setwd("/Users/keren/Desktop/covid_files/trees/")
tree_names = dir()
info = file.info(tree_names)
empty_files_plcs = which(info$size<500000)
tree_names = tree_names[-empty_files_plcs]
tmp = which(tree_names=="subs")
file_names = tree_names[-tmp]

init_seq = matrix(0,dim(regions)[1],1)
init_seq[all_no_subs_plcs] = (aligned_seqs[1][[1]][all_no_subs_plcs])
init_seq[which(init_seq==1)] = "A"
init_seq[which(init_seq==2)] = "C"
init_seq[which(init_seq==3)] = "G"
init_seq[which(init_seq==4)] = "T"

codon_names <- c("AAA","AAC","AAG","AAT","ACA","ACC","ACG","ACT","AGA","AGC","AGG","AGT","ATA","ATC","ATG",
                 "ATT","CAA","CAC","CAG","CAT","CCA","CCC","CCG","CCT","CGA","CGC","CGG","CGT","CTA","CTC",
                 "CTG","CTT","GAA","GAC","GAG","GAT","GCA","GCC","GCG","GCT","GGA","GGC","GGG","GGT","GTA",
                 "GTC","GTG","GTT","TAA","TAC","TAG","TAT","TCA","TCC","TCG","TCT","TGA","TGC","TGG","TGT",
                 "TTA","TTC","TTG","TTT")

get_node_names = function(tree_list){
  nodes_names = NULL
  for (i in 1:length(tree_list)){
    #print(i)
    nodes_names = rbind(nodes_names,cbind(tree_list[[i]]@name, tree_list[[i]]@edge))
  }
  colnames(nodes_names) = c("name","number")
  nodes_names = data.frame(nodes_names)
  nodes_names$number = as.numeric(nodes_names$number)
  return(nodes_names)
}
# load(tree_names[1])
# nodes_names = get_node_names(tree_list)
setwd("/Users/keren/Dropbox/covid/vars/")
# save(nodes_names,file = "nodes_names")
load(nodes_names)

get_seq = function(seq,tree_names,regions,missing_flag,curr_node){
  for (j in 1:length(tree_names)){
    load(tree_names[j])
    print(j)
    tmp = strsplit(tree_names[j],split = "")[[1]]
    site_pattern = as.numeric(paste0(tmp[11:length(tmp)],collapse = ""))
    if (length(which(tree_list[[curr_node]]@updated_probs==
                     max(tree_list[[curr_node]]@updated_probs)))>1){
      seq[which(regions$site_pattern==(site_pattern*missing_flag)),1] = "?"
    }else{
      seq[which(regions$site_pattern==(site_pattern*missing_flag)),1] = 
        tree_list[[curr_node]]@base
    }
  }
  return(seq)
}

get_node_details = function(parent_seq,parent_name,childe_name){
  # curr_edge = tree$edge[i,]
  # parent_name = nodes_names$name[which(nodes_names$number==curr_edge[1])]
  # child_name = nodes_names$name[which(nodes_names$number==curr_edge[2])]
  # edge_length = tree$edge.length[i]
  parent_node = nodes_names$number[which(nodes_names$name == parent_name)]
  child_node = nodes_names$number[which(nodes_names$name == parent_name)]
  setwd("/Users/keren/Desktop/covid_files/trees/")
  start_time = Sys.time()
  seq = get_seq(seq=init_seq,tree_names = tree_names,regions = regions,
                missing_flag = 1, curr_node = parent_name)
  end_time = Sys.time()
  setwd("/Users/keren/Desktop/covid_files/missing_trees/")
  start_time = Sys.time()
    seq = get_seq(seq = seq,tree_names = missing_tree_names,
                regions = regions,missing_flag = -1, curr_node = parent_name)
    end_time = Sys.time()
    print(end_time-start_time)
    L.N. = c("?",seq[1:(dim(regions)[1]-1)])
  R.N. = c(seq[2:dim(regions)[1]],"?")
  codon = matrix(0,dim(regions)[1],1)
  
}