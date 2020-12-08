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
# setwd("/Users/keren/Desktop/covid_files/trees/")
# load("tree_list_1")
# nodes_names = get_node_names(tree_list)
# setwd("/Users/keren/Dropbox/covid/vars/")
# save(nodes_names,file = "nodes_names")

setwd("/Users/keren/Dropbox/covid/vars/")
load("nodes_names")
load("site_details")
setwd("/Users/keren/Dropbox/EC2_files/")
source("COVID_functions.R")
tree = read.tree("tree.nwk")
setwd("/Users/keren/Dropbox/covid/vars/big_trees")
load("all_big_trees3")

#############
# remove incomplete codons:
codon_plcs = which(regions$codon.pos!=0)
plcs = which(diff(codon_plcs)!=1)
for (i in 1: length(plcs)){
  print(regions[codon_plcs[(plcs[i]-2):(plcs[i]+2)],])
}
gaps_in_coding = c(11118:11120,11123:11125,11132:11134,11137:11139,
                   21441:21443,22411:22419,25759:25761,26156:26158,
                   28324:28326,29646:29647)
plcs = NULL
for (i in 1:length(insertions)){
  if (sum(is.element(insertions[[i]],gaps_in_coding))>0)
    plcs = c(plcs,i)
}
plcs = NULL
for (i in 1:length(deletions)){
  print(i)
  if (sum(is.element(deletions[[i]],codon_plcs))>0)
    plcs = c(plcs,i)
}
############
# Functions:
update_states_and_outputs = function(curr_states,states,curr_outputs,
                                     outputs,uncertainty_plcs){
  if (length(uncertainty_plcs)>0){
    curr_states = curr_states[-uncertainty_plcs]
    curr_outputs = curr_outputs[-uncertainty_plcs,]
  }
  row.names(curr_outputs) = curr_states
  tmp = is.element(curr_states,states)
  new_states_plcs = which(tmp==FALSE)
  old_states_plcs = which(tmp==TRUE)
  if(length(old_states_plcs)>0){
    outputs[curr_states[old_states_plcs],] = 
      outputs[curr_states[old_states_plcs],] + curr_outputs[old_states_plcs,]
    # for (k in 1:length(old_states_plcs)){
    #   print(k)
    #   plc3 = which(states==curr_states[old_states_plcs[k]])
    #   outputs[plc3,] = outputs[plc3,] + curr_outputs[old_states_plcs[k],]
    # }
  }
  outputs = rbind(outputs,curr_outputs[new_states_plcs,])
  states = c(states,curr_states[new_states_plcs])
  return(list(states = states,outputs = outputs))
}
##########


# setwd("/Users/keren/Dropbox/mtDNA_tree_Build_17/")
# load("data.codon_nbs")
insertions = deletions =  list()
start_time = Sys.time()
for (i in 1:dim(tree$edge)[1]){
  print(i)
  curr_edge = tree$edge[i,]
  parent_name = nodes_names$name[which(nodes_names$number==curr_edge[1])]
  child_name = nodes_names$name[which(nodes_names$number==curr_edge[2])]
  # edge_length = tree$edge.length[i]
  parent_seq = strsplit(all_tree_list[[parent_name]]@probs,split = "")[[1]]
  child_seq = strsplit(all_tree_list[[child_name]]@probs,split = "")[[1]]
  insertion_plcs = which(parent_seq=="-" & (child_seq=="A" | 
                      child_seq=="C" | child_seq=="G" | child_seq=="T"))
  deletion_plcs = which(child_seq=="-" & (parent_seq=="A" | 
                        parent_seq=="C" | parent_seq=="G" | parent_seq=="T"))
  is.element(insertion_plcs,codon_plcs)
  deletions[[i]] = deletion_plcs
  seq_inserted
}
end_time = Sys.time()
print(end_time-start_time)

base_outputs = base_outputs[-1,]
base_states = base_states[-1]
codon_states = codon_states[-1]
codon_outputs = codon_outputs[-1,]

setwd("/Users/keren/Dropbox/covid/vars/")
save(base_states,file = "base_states")
save(codon_states,file = "codon_states")
save(base_outputs,file = "base_outputs")
save(codon_outputs,file = "codon_outputs")
