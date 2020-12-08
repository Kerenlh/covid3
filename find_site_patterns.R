find_site_patterns = function(parallel_nums){
  # setwd("/Users/keren/Dropbox/covid")
  # setwd("C:/Users/ifog/Dropbox/covid")
  # setwd("/home/ubuntu/Dropbox/covid")
  # setwd("/Users/keren/Dropbox/covid/vars")
  # setwd("C:/Users/ifog/Dropbox/covid/vars")
  # setwd("/home/ubuntu/Dropbox/EC2_files")
  setwd("/Users/keren/Dropbox/EC2_files")
  source("COVID_functions.R")
  all_no_subs_plcs = read.csv("all_no_subs_plcs.csv")
  all_no_subs_plcs = as.numeric(all_no_subs_plcs[,2])
  aligned_seqs = readDNAStringSet(file="aligned_tree_seqs_decipher")
  tree = read.tree("tree.nwk")
  tree = short_names_tree(tree)
  
  site_patterns = NULL
  no_subs_plcs = NULL
  site_patterns_plcs = NULL
  # setwd("/Users/keren/Dropbox/EC2_files/trees/")
  # setwd("/home/ubuntu/Dropbox/EC2_files/trees/")
  setwd("/Users/keren/Dropbox/covid/vars/missing_site_patterns2/")
  missing_plcs = c(9857, 9858, 9859, 9862, 9863, 9864, 9865, 9866, 9869,
                   9871, 9873, 9879, 9881, 9891, 9899, 9900)
  # for (i in all_no_subs_plcs[]){
  # for (i in all_no_subs_plcs[((parallel_nums-1)*10+1) : min((parallel_nums*10), length(aligned_seqs[[1]]))  ]){ #(i in :length(aligned_seqs[[1]])){
  for (i in missing_plcs){  
    print(i)
    curr_site_pattern = NULL
    start_time = Sys.time()
    for (j in 1:length(aligned_seqs)){
      curr_site_pattern = c(curr_site_pattern, as.character(aligned_seqs[j][[1]][i]))
    }
    end_time = Sys.time()
    print(end_time-start_time)
    if ( length(table(curr_site_pattern))==1) { 
      # |(length(table(curr_site_pattern))==2 & is.element("-",curr_site_pattern))){
      no_subs_plcs = c(no_subs_plcs,i)
    }else{
      site_patterns = c(site_patterns,paste0(curr_site_pattern, collapse = ""))
      site_patterns_plcs = c(site_patterns_plcs,i) 
    }
    save(site_patterns,file = paste0("site_patterns_",i,collapse = ""))
    save(site_patterns_plcs,file = paste0("site_patterns_plcs_",i,collapse = ""))
    save(no_subs_plcs,file = paste0("no_subs_plcs_",i,collapse = ""))
  }
  # setwd("/home/ubuntu/Dropbox/EC2_files/trees/")
  return(site_patterns)
}

# save(site_patterns, file = "site_patterns")
# load("site_patterns")

# all_site_patterns = unique(site_patterns)
# site_pattern_indices = matrix(0,length(site_patterns),1)
# for (i in 1:length(all_site_patterns)){
#   site_pattern_indices[which(site_patterns==all_site_patterns[i])] = i
# }
require(snow)
library(snow)
library(doSNOW)
no_cores = 64
cl = makeCluster(no_cores,type="SOCK",outfile="")
registerDoSNOW(cl)
parallel_nums = c(1:819)
print(Sys.time())
site_patterns = clusterApply(cl, parallel_nums,find_site_patterns)
print(Sys.time())
stopCluster(cl)
