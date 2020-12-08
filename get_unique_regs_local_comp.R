# setwd("/a/home/cc/math/kerenlh/EC2_files/Runs/unique_regs/i_3_9")
# wait = runif(1, 0, 30)
# Sys.sleep(wait)
# if (file.exists("begin")){
#   load("begin")
# }else{
#   begin = 1
# }
# 
# old_begin = begin
# 
# job_length = 500#10000 #last-begin+1
# last = begin+job_length-1
# begin = last+1
# save(begin,file = "begin")
# 
# begin = old_begin
# print(begin)
# print(last)

# setwd("/Users/keren/Dropbox/covid/vars/datasets/")
# load("missing_regs_dataset2_i_3_9")
# dataset = missing_regs_dataset

setwd("/Users/keren/Dropbox/covid/vars/datasets/syn_i_1/")
load("missing_dataset_i_1")
dataset = missing_dataset
plcs = which(duplicated(dataset$reg_name2)==FALSE)
dataset = dataset[plcs,]
dataset$output = 1

# setwd("/a/home/cc/math/kerenlh/EC2_files")
setwd("/Users/keren/Dropbox/covid/vars/")
load("codons_table3")
load("iterate_vals_codons")
source("/Users/keren/Dropbox/covid/debug_functions.R")
# source("/a/home/cc/math/kerenlh/EC2_files/debug_functions.R")
library(hash)
library(digest)

reg_hash = hash()

for (i in 1:dim(dataset)[1]){
  start_time = Sys.time()
  print(i)
  tmp = get_data(dataset[i,],iterate_vals = iterate_vals)
  # curr_data = tmp[[1]]; 
  curr_details = tmp[[2]]
  reg_hash[[dataset$reg_name2[i]]] = 
    curr_details[(length(iterate_vals)+3):length(curr_details)]
  print(Sys.time()-start_time)
}

# setwd("/a/home/cc/math/kerenlh/EC2_files/Runs/unique_regs/i_3_9")
setwd("/Users/keren/Dropbox/covid/vars/datasets/reg_hashes/i_1/")
save(reg_hash,file = paste0("missing_reg_hash_i_1",collapse = ""))
stop("so that the workspace is not saved")

