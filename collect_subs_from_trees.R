collect_subs = function(parallel_nums){
  setClass(Class='Tree', 
           representation(name= 'character',edge = 'numeric', parent_edge = 'numeric',
                          parent_name = 'character',  
                          children_edges = 'matrix', children_names = 'character', 
                          probs = 'matrix', updated_probs = 'matrix', base = 'character',leaf = 'logical'))
  
  collect_subs_from_tree = function(tree_list,curr_name,subs,before_edge,u_flag){
    #print(curr_name)
    curr_edge = tree_list[[curr_name]]@edge
    parent_name = tree_list[[curr_name]]@parent_name
    after = tree_list[[curr_name]]@base
    if (length(which(tree_list[[curr_name]]@updated_probs==
                     max(tree_list[[curr_name]]@updated_probs)))>1){
      u_flag = 10*u_flag+1 
    }
    before = tree_list[[parent_name]]@base
    if (before!=after){
      new_sub = c(before,after,before_edge,curr_edge,u_flag,tree_list[[curr_name]]@leaf)
      subs = rbind(subs,new_sub)
      print(new_sub)
      curr_edges = NULL
      u_flag = 0
      before_edge = curr_edge
    }
    curr_children = tree_list[[curr_name]]@children_names
    if (tree_list[[curr_name]]@leaf==FALSE){
      for (i in 1:length(curr_children)){
        subs = collect_subs_from_tree(tree_list = tree_list,
                                      curr_name = curr_children[i], subs = subs,
                                      before_edge = before_edge, u_flag = u_flag)
      }
    }
    return(subs)
  }
  setwd("/Users/keren/Dropbox/covid/vars/")
  load("missing_subs")
  root_edge = 12225
  root_name = paste0(c("no_name_",root_edge),collapse = "")
  # Define parent for root:
  tt = new("Tree", name = "no_name_0", edge = 0, parent_edge = 0,
           parent_name = "root_grandpa", children_edges = matrix(12225,1,1),
           children_names = "no_name_12225", base = "0", probs = matrix(0,1,5),
           updated_probs = matrix(0,1,5), leaf = FALSE)
  setwd("/Users/keren/Desktop/covid_files/trees")
  # setwd("/home/ubuntu/Dropbox/EC2_files")
  file_names = dir()
  info = file.info(file_names)
  empty_files_plcs = which(info$size<500000)
  file_names = file_names[-empty_files_plcs]
  # tmp = which(file_names=="subs")
  # file_names = file_names[-tmp]
  for (i in missing_subs[parallel_nums]){#1:length(file_names)){
    curr_file_name = paste0("tree_list_",i)
    load(curr_file_name)
    if (length(tree_list)==0){
      next
    }
    print(curr_file_name)
    # tmp = strsplit(file_names[i],split = "")[[1]]
    # site_pattern = as.numeric(paste0(tmp[11:length(tmp)],collapse=""))
    site_pattern = i
    subs = NULL
    u_flag = 0
    curr_edges = NULL
    tree_list[["no_name_0"]]=tt
    start_time = Sys.time()
    subs = collect_subs_from_tree(tree_list = tree_list,curr_name = root_name,
                                  subs = subs,before_edge = root_edge,u_flag = u_flag)
    print(Sys.time()-start_time)
    colnames(subs) = c("before","after","before_edge","after_edge","uncertainty_flag","leaf")
    new_subs = cbind(subs,site_pattern)
    new_subs = new_subs[-1,]
    # load("subs")
    subs = new_subs
    setwd("/Users/keren/Desktop/covid_files/subs/")
    save(subs,file = paste0("subs_",site_pattern,collapse = ""))
  }
}
require(snow)
library(snow)
library(doSNOW)
no_cores = 2
cl = makeCluster(no_cores,type="SOCK",outfile="")
registerDoSNOW(cl)
#cluster_delta = 50 # should update inside the function as well
#length(unique_site_patterns)=13366
parallel_nums = c(1:1103)
print(Sys.time())
setwd("/Users/keren/Dropbox/covid/vars")
# tree = read.tree("tree.nwk")
# tree = short_names_tree(tree)
site_patterns = clusterApply(cl, parallel_nums,collect_subs)
print(Sys.time())
stopCluster(cl)

