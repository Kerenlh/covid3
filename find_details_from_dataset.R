library(hash)
library(digest)

get_all_reg_hash = function(file_names){
  all_reg_hash = hash()
  for (i in 1:length(file_names)){
    if (paste0(strsplit(file_names[i],split = "")[[1]][1:9]
               ,collapse = "")=="reg_hash."){ # change to reg_hash.
      load(file_names[i])
      curr_hash_names = names(reg_hash)
      print(i)
      for (j in 1:length(curr_hash_names)){
        all_reg_hash[[curr_hash_names[j]]] = reg_hash[[curr_hash_names[j]]]
      }
    }
  }
  return(all_reg_hash)
}
setwd("/Users/keren/Dropbox/covid/vars/datasets/reg_hashes/base/")
file_names = dir()
all_reg_hash = get_all_reg_hash(file_names = file_names)
save(all_reg_hash,file = "all_reg_hashes")

setwd("/Users/keren/Dropbox/covid/vars/datasets/base/reg_hash_plcs/")
load("all_reg_hash_plcs_base")

setwd("/Users/keren/Dropbox/covid/vars/datasets/reg_hashes/base/")
load("all_reg_hashes")
setwd("/Users/keren/Dropbox/covid/vars/datasets/")
load("dataset_reg_name_base")
length_reg_hash = 18
reg_output = matrix(0,length(dataset_reg_name),length_reg_hash)
reg_hash_names = names(all_reg_hash)
print(length(reg_hash_names))
for (j in 1:length(reg_hash_names)){
  print(j)
  plcs = all_reg_hash_plcs[[j]]
  if (length(plcs)>0){
    reg_output[plcs,] = matrix(all_reg_hash[[reg_hash_names[j]]],
                               length(plcs),dim(reg_output)[2],byrow = TRUE)
  }
}
setwd("/Users/keren/Dropbox/covid/vars/datasets/reg_hashes/base/")
save(reg_output,file = "reg_output_base")
rm(reg_output)
outputs_list = list()
outputs_list[[1]] = c(1,2) #syn,non_syn
outputs_list[[2]] = c(3,4) # transition,transversion
outputs_list[[3]] = 5 # y
outputs_list[[4]] = c(6:9) # A/C/G/T
for (i in 2:length(outputs_list)){
  setwd("/Users/keren/Dropbox/covid/vars/datasets/")
  load("dataset_base")
  dataset = dataset_base
  rm(dataset_base)
  plcs = which(is.element(as.numeric(dataset$output),outputs_list[[i]])==TRUE)
  details = dataset[,1:(dim(dataset)[2]-1)]
  rm(dataset)
  setwd("/Users/keren/Dropbox/covid/vars/datasets/reg_hashes/base/")
  load("reg_output_base")
  details = cbind(details,matrix(0,dim(details)[1],1),reg_output[plcs,])
  colnames(details) = c("L.neighbor","R.neighbor","codon",            
                        "codon.pos","gene","mat_peptide","stem_loop","CG","amino_acid",       
                        "output","ID","logLik","df","AIC","theta",
                        "NB.converged","rows_num","cols_num","P.GLR","P.test","P.logLik",
                        "P.df","P.AIC","P.converged","num_non_0_rows_1","output_sum",
                        "added_log_lik","added_log_lik.P", "output_0_or_1_row") 
  setwd("/Users/keren/Dropbox/covid/details//")
  save(details,file = paste0("details_base.",i,collapse = ""))
}

