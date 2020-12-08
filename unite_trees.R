unite_trees = function(parallel_nums){
  cluster_delta = 1000
  print(paste("parallel_nums = ",parallel_nums))
  setwd("/Users/keren/Desktop/covid_files/missing_trees/")
  missing_tree_names = dir()
  # setwd("/Users/keren/Desktop/covid_files/trees/")
  # tree_names = dir()
  # info = file.info(tree_names)
  # empty_files_plcs = which(info$size<500000)
  # tree_names = tree_names[-empty_files_plcs]
  # tmp = which(tree_names=="subs")
  # file_names = tree_names[-tmp]
  
  setwd("/Users/keren/Dropbox/covid/vars/")
  load("missing_all_no_subs_plcs")
  # setwd("/Users/keren/Dropbox/covid/vars/missing_site_patterns2/trees/")
  # missing_tree_names = dir()
  
  setClass(Class='Tree',
           representation(name= 'character',edge = 'numeric', parent_edge = 'numeric',
                          parent_name = 'character',
                          children_edges = 'matrix', children_names = 'character',
                          probs = 'character', updated_probs = 'matrix', base = 'character',leaf = 'logical'))
  # tt = new("Tree", name = "a", edge = 1, parent_edge = 2,
  #          parent_name = "b", children_edges = matrix(1,1,2),
  #          children_names = "c", base = "t", probs = matrix(0,1,2),
  #          updated_probs = matrix(0,1,2), leaf = TRUE)
  # tree_list[[curr_name]] = tt
  
  setwd("/Users/keren/Dropbox/EC2_files/")
  # setwd("/home/ubuntu/Dropbox/EC2_files")
  load("site_details")
  source("COVID_functions.R")
  # load("new_tree_list")
  aligned_seqs = readDNAStringSet(file="aligned_tree_seqs_decipher")
  init_seq = matrix("M",length(aligned_seqs[1][[1]]),1)
  init_seq[all_no_subs_plcs] = (aligned_seqs[1][[1]][all_no_subs_plcs])
  init_seq[which(init_seq==1)] = "A"
  init_seq[which(init_seq==2)] = "C"
  init_seq[which(init_seq==3)] = "G"
  init_seq[which(init_seq==4)] = "T"
  init_seq = paste0(init_seq,collapse = "")
  
  
  setwd("/Users/keren/Desktop/covid_files/missing_trees/")
  # setwd("/Users/keren/Desktop/covid_files/trees/")
  load("tree_list_1")
  new_tree_list = tree_list
  for (i in 1:length(new_tree_list)){
    # print(i)
    new_tree_list[[i]]@probs = init_seq
  }
  # save(new_tree_list,file = "new_tree_list")
  
  # setwd("/Users/keren/Dropbox/covid/vars/missing_site_patterns2/trees/")
  # setwd("/Users/keren/Dropbox/EC2_files/trees/")
  # setwd("/home/ubuntu/Dropbox/EC2_files/trees/")
  # setwd("/Users/keren/Desktop/covid_files/missing_trees/")
  tree_names = missing_tree_names
  # tmp = which(tree_names=="subs" | tree_names=="nodes_names")
  # tree_names = tree_names[-tmp]
  missing_flag = -1
  start1 = Sys.time()
  # for (i in c((1+(parallel_nums-1)*cluster_delta):(parallel_nums*cluster_delta))){
    for (i in 1:length(tree_names)){
    print(i)
    load(tree_names[i])
    tree_list[["no_name_12225"]]@updated_probs = tree_list[["no_name_12225"]]@probs
    tmp = strsplit(tree_names[i],split = "")[[1]]
    site_pattern = as.numeric(paste0(tmp[11:length(tmp)],collapse = ""))
    plcs = which(regions$site_pattern==(site_pattern*missing_flag))
    start_time = Sys.time()
    for(j in 1:length(new_tree_list)){
      if (length(which(tree_list[[j]]@updated_probs==
                       max(tree_list[[j]]@updated_probs)))>1){
        # tmp = strsplit(new_tree_list[[j]]@probs,split = "")[[1]]
        # tmp[plcs] = "?"
        # new_tree_list[[j]]@probs = paste0(tmp,collapse = "")
        for (k in 1:length(plcs)){
          substr(new_tree_list[[j]]@probs,plcs[k],plcs[k]) = "?"
          # print(paste(tree_names[i],"j=",j))
        }
      }else{
        # tmp = strsplit(new_tree_list[[j]]@probs,split = "")[[1]]
        # tmp[plcs] = tree_list[[j]]@base
        # new_tree_list[[j]]@probs = paste0(tmp)
        # 
        for (k in 1:length(plcs)){
          substr(new_tree_list[[j]]@probs,plcs[k],plcs[k]) =  tree_list[[j]]@base
        }
      }
    }
    end_time = Sys.time()
    # print(end_time-start_time)
  }
  print(paste("1000 runing time = ",Sys.time()-start1))
  setwd("/Users/keren/Dropbox/covid/vars/big_trees")
  # setwd("/home/ubuntu/Dropbox/EC2_files")
  # save(new_tree_list,file = paste0(c("new_tree_list_",parallel_nums),collapse = ""))
  save(new_tree_list,file = paste0(c("new_tree_list_",15),collapse = ""))
}
# 
require(snow)
library(snow)
library(doSNOW)
no_cores = 3
cl = makeCluster(no_cores,type="SOCK",outfile="")
registerDoSNOW(cl)
cluster_delta = 1000 # should update inside the function as well
# length(unique_site_patterns)=13366
setwd("/Users/keren/Dropbox/EC2_files/trees/")
# setwd("/home/ubuntu/Dropbox/EC2_files/trees/")
tree_names = dir()
parallel_nums = c(1:14)#(13366/cluster_delta))
# print(Sys.time())
site_patterns = clusterApply(cl, parallel_nums,unite_trees)
# print(Sys.time())
stopCluster(cl)


