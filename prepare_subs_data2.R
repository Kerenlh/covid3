prepare_data = function(parallel_nums){
  cluster_delta = 5027
  
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
  gaps_in_coding = c(11118:11120,11123:11125,11132:11134,11137:11139,
                     21441:21443,22411:22419,25759:25761,26156:26158,
                     28324:28326)
  # 29646:29647: This is also a gap of 2, an insertion appears only in one seq, 
  # relevant node is i=5350,  parent_name = "no_name_14339", child_name "MT520337.1"
  
  regions$codon.pos[gaps_in_coding] = matrix(c(1:3),33,1)
  
  regions$codon.pos[11123:11125] = c(3,1,2)
  regions$codon.pos[11132:11134] = c(3,1,2)
  regions$codon.pos[11137:11139] = c(2,3,1)
  regions$codon.pos[22411:22419] = c(2,3,1,2,3,1,2,3,1)
  regions$codon.pos[21612:21613] = 0
  regions$codon.pos[27827:27828] = 0
  regions$codon.pos[25759:25761] = c(3,1,2)
  
  codon_pos_1_plcs = which(regions$codon.pos==1)
  codon_pos_2_plcs = which(regions$codon.pos==2)
  codon_pos_3_plcs = which(regions$codon.pos==3)
  codon_pos_plcs = cbind(codon_pos_1_plcs,codon_pos_2_plcs,codon_pos_3_plcs)
  colnames(codon_pos_plcs) = c("pos1","pos2","pos3")
  codon_pos_plcs = data.frame(codon_pos_plcs)
  tmp = which(codon_pos_2_plcs-codon_pos_1_plcs!=1)
  tmp = which(codon_pos_3_plcs-codon_pos_1_plcs!=2)
  
  codon_plcs = which(regions$codon.pos!=0)
  
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
  
  base_states = "tmp"
  codon_states = "tmp"
  init_outputs = matrix(0,dim(regions)[1],10)
  colnames(init_outputs) = c("A","C","G","T","line","q","s","M","exposure","zero_branch")
  init_outputs = data.frame(init_outputs)
  base_outputs = init_outputs[1,]
  codon_outputs = init_outputs[1,]
  edges_with_no_subs = NULL
  
  # for (i in c((1+(parallel_nums-1)*cluster_delta):(parallel_nums*cluster_delta))){
  for (i in 15082:dim(tree$edge)[1]){
    # start_time = Sys.time()
    if (i==5350){
      next
    }
    print(i)
    curr_edge = tree$edge[i,]
    parent_name = nodes_names$name[which(nodes_names$number==curr_edge[1])]
    child_name = nodes_names$name[which(nodes_names$number==curr_edge[2])]
    edge_length = tree$edge.length[i]
    parent_seq = strsplit(all_tree_list[[parent_name]]@probs,split = "")[[1]]
    child_seq = strsplit(all_tree_list[[child_name]]@probs,split = "")[[1]]
    subs_plcs = which((parent_seq==child_seq)==FALSE)
    curr_outputs = init_outputs
    if(length(subs_plcs)>0){
      for (j in 1:length(subs_plcs)){
        if (child_seq[subs_plcs[j]]=="-"){
          curr_outputs[subs_plcs[j],"line"]=1
        }else if (child_seq[subs_plcs[j]]=="?"){
          curr_outputs[subs_plcs[j],"q"]=1
        }else if (child_seq[subs_plcs[j]]=="*"){
          curr_outputs[subs_plcs[j],"s"]=1
        }else{
          curr_outputs[subs_plcs[j],child_seq[subs_plcs[j]]] = 1
        }
      }
    }else{
      edges_with_no_subs = c(edges_with_no_subs,i)
    }
    curr_outputs$exposure = edge_length
    if (edge_length==0){
      curr_outputs$zero_branch = 1
    }
    
    gaps_plcs = which(parent_seq=="-")
    short_parent = parent_seq[-gaps_plcs]
    curr_codon_plcs = codon_plcs[-which(is.element(codon_plcs,gaps_plcs))]
    codon_short_parent = parent_seq[curr_codon_plcs]
    curr_sites = regions$site[-gaps_plcs]
    curr_codon_sites = regions$site[curr_codon_plcs]
    curr_codon.pos = regions$codon.pos[curr_codon_plcs]
    codons = NULL
    k=1
    while(k <=length(curr_codon_sites)){
      if (sum(curr_codon.pos[k:(k+2)]==c(1,2,3))==3){
        codons = c(codons,rep(paste0(codon_short_parent[k:(k+2)],collapse = ""),3))
        k = k+3
      }else{
        codons = c(codons,"missing")
        print(paste("k=",k,"i=",i))
        k= k+1
      }
    }
    L.N. = c("*",short_parent[1:(length(short_parent)-1)])
    R.N. = c(short_parent[2:length(short_parent)],"*")
    L.N.codon = parent_seq[curr_codon_plcs-1]
    R.N.codon = parent_seq[curr_codon_plcs+1]
    # codons = paste0(short_parent[codon_pos_1_plcs],
    #                 short_parent[codon_pos_2_plcs],
    #                 short_parent[codon_pos_3_plcs])
    # codons = matrix(codons,length(codons),3)
    # codons = matrix(t(codons),dim(codons)[1]*3,1,byrow = TRUE)
    curr_codon_states = paste0(L.N.codon,R.N.codon,
                               parent_seq[curr_codon_plcs],codons,
                               regions$site[curr_codon_plcs])
    codon_uncertainty_plcs = which(parent_seq[curr_codon_plcs]=="?" | 
                                     curr_outputs$q[curr_codon_plcs]=="1")
    tmp = update_states_and_outputs(curr_states = curr_codon_states,
                                    states = codon_states,
                                    curr_outputs = curr_outputs[curr_codon_plcs,],
                                    outputs = codon_outputs,
                                    uncertainty_plcs = codon_uncertainty_plcs)
    codon_states = tmp$states
    codon_outputs = tmp$outputs
    
    curr_base_states = paste0(L.N.,R.N.,short_parent,curr_sites)
    base_uncertainty_plcs = which(parent_seq=="?" | 
                                    curr_outputs$q=="1")
    tmp = update_states_and_outputs(curr_states = curr_base_states,
                                    states = base_states,
                                    curr_outputs = curr_outputs[curr_sites,],
                                    outputs = base_outputs,
                                    uncertainty_plcs = base_uncertainty_plcs)
    base_states = tmp$states
    base_outputs = tmp$outputs
    # end_time = Sys.time()
    # print(end_time-start_time)
    
  }
  
  base_outputs = base_outputs[-1,]
  base_states = base_states[-1]
  codon_states = codon_states[-1]
  codon_outputs = codon_outputs[-1,]
  
  setwd("/Users/keren/Dropbox/covid/vars/")
  save(base_states,file = paste0("base_states_",parallel_nums,collapse =""))
  save(codon_states,file = paste0("codon_states_",parallel_nums,collapse =""))
  save(base_outputs,file =  paste0("base_outputs_",parallel_nums,collapse =""))
  save(codon_outputs,file = paste0("codon_outputs_",parallel_nums,collapse =""))
}
require(snow)
library(snow)
library(doSNOW)
no_cores = 1
cl = makeCluster(no_cores,type="SOCK",outfile="")
registerDoSNOW(cl)
# codon_vals = 37:66
parallel_nums = c(4)
print(Sys.time())
details = clusterApply(cl, parallel_nums,prepare_data)
print(Sys.time())
stopCluster(cl)

