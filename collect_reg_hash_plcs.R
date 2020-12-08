setwd("/Users/keren/Dropbox/covid/vars/datasets/base/reg_hash_plcs/")
file_names = dir()
all_reg_hash_plcs = list()
for (i in 1:length(file_names)){
  load(file_names[i])
  print(i)
  for(j in 1:length(reg_hashes_plcs)){
    if (is.null(reg_hashes_plcs[[j]])==FALSE){
      # print(j)
      all_reg_hash_plcs[[j]] = reg_hashes_plcs[[j]]
    }
  }
}
save(all_reg_hash_plcs,file = "all_reg_hash_plcs_base")