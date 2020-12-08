setwd("/a/home/cc/math/kerenlh/EC2_files/Runs/datasets/reg_hashes/base")
# setwd("/Users/keren/Dropbox/EC2_files/")
wait = runif(1, 0, 30)
Sys.sleep(wait)
if (file.exists("reg_hash_num")){
  load("reg_hash_num")
}else{
  reg_hash_num = 1
}

old_reg_hash_num = reg_hash_num
iter_num = 1000
reg_hash_num = reg_hash_num+iter_num
save(reg_hash_num,file = "reg_hash_num")

reg_hash_num = old_reg_hash_num
print(reg_hash_num)


# setwd("/Users/keren/Dropbox/covid/vars/datasets/reg_hashes/base/")
load("all_reg_hashes")
# setwd("/Users/keren/Dropbox/covid/vars/datasets/")
# load("dataset_base")
# dataset_reg_name = dataset_base$reg_name2
# save(dataset_reg_name, file = "dataset_reg_name_base")
# rm(dataset_base)
library(hash)
load("dataset_reg_name_base")
reg_hash_names = names(all_reg_hash)
print(length(reg_hash_names))
reg_hashes_plcs = hash()
for (j in reg_hash_num:(reg_hash_num+iter_num-1)){
  print(j)
  plcs = which(dataset_reg_name==reg_hash_names[j])
  reg_hashes_plcs[[ reg_hash_names[j] ]] = plcs
}
save(reg_hashes_plcs,file = paste0("reg_hashes_plcs_base.",reg_hash_num,collapse = ""))
stop("so that the workspace is not saved")


