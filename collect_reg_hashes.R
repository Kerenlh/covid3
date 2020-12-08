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

setwd("/Users/keren/Dropbox/covid/vars/datasets/reg_hashes/i_3_9/")
file_names = dir()
all_reg_hash = get_all_reg_hash(file_names = file_names)
# The hash is split to 2 because of memory issues
all_reg_hash = get_all_reg_hash(file_names = file_names[1:100])
all_reg_hash = get_all_reg_hash(file_names = file_names[101:199])
all_reg_hash = get_all_reg_hash(file_names = file_names[200:290])
all_reg_hash = get_all_reg_hash(file_names = file_names[291:length(file_names)])
save(all_reg_hash,file = "all_reg_hash1")
save(all_reg_hash,file = "all_reg_hash2")
save(all_reg_hash,file = "all_reg_hash3")
save(all_reg_hash,file = "all_reg_hash4")

setwd("/Users/keren/Dropbox/covid/vars/datasets/reg_hashes/i_3_9/missing_reg_hashes/")
file_names = dir()
all_reg_hash = get_all_reg_hash(file_names = file_names)
save(all_reg_hash,file = "missing_reg_hashes")

setwd("/Users/keren/Dropbox/covid/vars/datasets/reg_hashes/base/")
file_names = dir()
all_reg_hash = get_all_reg_hash(file_names = file_names)
save(all_reg_hash,file = "all_reg_hashes")

##### collect and save all reg_hashes to 1 hash


setwd("/Users/keren/Dropbox/covid/vars/datasets/i_3_9/")
file_names = dir()
dataset_3_4 = dataset_5 = dataset_6_9 = NULL
load(file_names[3000])
plc1 = (which(colnames(short_details)=="reg_name2"))
for (i in 1:length(file_names)){
  print(i)
  start_time = Sys.time()
  file_name_split = strsplit(file_names[i],split = "")[[1]]
  if (paste0(file_name_split[1:14]
             ,collapse = "")=="short_details."){
    load(file_names[i])
    dataset_3_4 = 
      rbind(dataset_3_4,short_details[which(short_details$output==3
                                            |short_details$output==4),1:plc1])
    dataset_5 = 
      rbind(dataset_5,short_details[which(short_details$output==5),1:plc1])
    dataset_6_9 = 
      rbind(dataset_6_9,short_details[which(short_details$output>=6 & 
                                              short_details$output<=9),1:plc1])
  }
  print(Sys.time()-start_time)
}
setwd("/Users/keren/Dropbox/covid/vars/datasets/")
save(dataset_3_4,file = "dataset_3_4")
save(dataset_5,file = "dataset_5")
save(dataset_6_9,file = "dataset_6_9")


setwd("/Users/keren/Dropbox/covid/vars/datasets/base//")
file_names = dir()
dataset_base = NULL
load(file_names[1])
plc1 = (which(colnames(short_details)=="reg_name2"))
for (i in 1:length(file_names)){
  print(i)
  start_time = Sys.time()
  file_name_split = strsplit(file_names[i],split = "")[[1]]
  if (paste0(file_name_split[1:14]
             ,collapse = "")=="short_details."){
    load(file_names[i])
    dataset_base = rbind(dataset_base,short_details[,1:plc1])
   }
  print(Sys.time()-start_time)
}
setwd("/Users/keren/Dropbox/covid/vars/datasets/")
save(dataset_base,file = "dataset_base")



####### Collect all dataset to 3 datasets divided according to the output


add_reg_hash_details_to_dataset = function(dataset,details_file_name){
  setwd("/Users/keren/Dropbox/covid/vars/datasets/reg_hashes/i_3_9/")
  details = NULL
  length_reg_hash = 15
  reg_output = matrix(0,dim(dataset)[1],length_reg_hash)
  for (i in 1:4){
    start_time = Sys.time()
    load(paste0("all_reg_hash",i,collapse = ""))
    reg_hash_names = names(all_reg_hash)
    print(length(reg_hash_names))
    print(i)
    for (j in 1:length(reg_hash_names)){
      # print(j)
      plcs = which(dataset$reg_name2==reg_hash_names[j])
      if (length(plcs)>0){
        reg_output[plcs,] = matrix(all_reg_hash[[reg_hash_names[j]]],
                                   length(plcs),dim(reg_output)[2],byrow = TRUE)
      }
    }
    print(paste("i=",i))
    print(Sys.time()-start_time)
  }
  details = cbind(dataset[,1:(dim(dataset)[2]-1)],matrix(0,dim(dataset)[1],1),reg_output)
  colnames(details) = c("L.neighbor","R.neighbor","codon",            
                        "codon.pos","gene","mat_peptide","stem_loop","CG","amino_acid",       
                        "output","ID","logLik","df","AIC","theta",
                        "NB.converged","rows_num","cols_num","P.GLR","P.test","P.logLik",
                        "P.df","P.AIC","P.converged","num_non_0_rows_1","output_sum") 
  setwd("/Users/keren/Dropbox/covid/details/details_i_3_9/")
  save(details,file = details_file_name)
  return(details)
}

setwd("/Users/keren/Dropbox/covid/vars/datasets/")

load("dataset_3_4")
details_3_4 = add_reg_hash_details_to_dataset(dataset = dataset_3_4,
                                              details_file_name = "details_3_4")

load("dataset_5")
details_5 = add_reg_hash_details_to_dataset(dataset = dataset_5,
                                            details_file_name = "details_5")

load("dataset_6_9")
details_6_9 = add_reg_hash_details_to_dataset(dataset = dataset_6_9,
                                              details_file_name = "details_6_9")

#########################
# Get missing regs dataset: 
get_missing_regs_dataset = function(details,dataset){
  plcs = which(details$rows_num==0)
  missing_regs_dataset = dataset[plcs,]
  plcs2 = which(duplicated(missing_regs_dataset$reg_name2)==FALSE)
  missing_regs_dataset = missing_regs_dataset[plcs2,]
  rm(dataset)
  rm(details)
  return(missing_regs_dataset)
}

setwd("/Users/keren/Dropbox/covid/details/details_i_3_9/")
load("details_3_4")
setwd("/Users/keren/Dropbox/covid/vars/datasets/")
load("dataset_3_4")
missing_regs_dataset_3_4 = get_missing_regs_dataset(details,dataset_3_4)
rm(dataset_3_4)
rm(details)
setwd("/Users/keren/Dropbox/covid/details/details_i_3_9/")
load("details_5")
setwd("/Users/keren/Dropbox/covid/vars/datasets/")
load("dataset_5")
missing_regs_dataset_5 = get_missing_regs_dataset(details,dataset_5)
rm(dataset_5)
rm(details)
setwd("/Users/keren/Dropbox/covid/details/details_i_3_9/")
load("details_6_9")
setwd("/Users/keren/Dropbox/covid/vars/datasets/")
load("dataset_6_9")
missing_regs_dataset_6_9 = get_missing_regs_dataset(details,dataset_6_9)
rm(dataset_6_9)
rm(details)

missing_regs_dataset = rbind(missing_regs_dataset_3_4,
                             missing_regs_dataset_5,
                             missing_regs_dataset_6_9)
plcs = which(duplicated(missing_regs_dataset$reg_name2)==FALSE)
missing_regs_dataset = missing_regs_dataset[plcs,]
setwd("/Users/keren/Dropbox/covid/vars/datasets/")
save(missing_regs_dataset,file = "missing_regs_dataset2_i_3_9")


#####
# Add reg hashes details to the missing dataset:
setwd("/Users/keren/Dropbox/covid/vars/datasets/reg_hashes/i_3_9/missing_reg_hashes/")
load("missing_reg_hashes")
setwd("/Users/keren/Dropbox/covid/details/details_i_3_9/")
load("details_3_4_missing_likelihood")
setwd("/Users/keren/Dropbox/covid/vars/datasets/")
load("dataset_3_4")
dataset = dataset_3_4
rm(dataset_3_4)
tmp = apply(dataset[,1:10],1,paste0,collapse = ".")
tmp2 = apply(details[,1:10],1,paste0,collapse = ".")
missing_plcs = which(is.element(tmp,tmp2)==FALSE)
dim(dataset)[1]-dim(details)[1]
missing_dataset = dataset[missing_plcs,]

dataset = missing_dataset
length_reg_hash = 15
reg_output = matrix(0,dim(dataset)[1],length_reg_hash)
reg_hash_names = names(all_reg_hash)
print(length(reg_hash_names))
for (j in 1:length(reg_hash_names)){
  print(j)
  plcs = which(dataset$reg_name2==reg_hash_names[j])
  if (length(plcs)>0){
    reg_output[plcs,] = matrix(all_reg_hash[[reg_hash_names[j]]],
                               length(plcs),dim(reg_output)[2],byrow = TRUE)
  }
}
setwd("/Users/keren/Dropbox/covid/vars/datasets/reg_hashes/")
load("missing_reg_hash2")
all_reg_hash = reg_hash
reg_hash_names = names(all_reg_hash)
print(length(reg_hash_names))
for (j in 1:length(reg_hash_names)){
  print(j)
  plcs = which(dataset$reg_name2==reg_hash_names[j])
  if (length(plcs)>0){
    reg_output[plcs,] = matrix(all_reg_hash[[reg_hash_names[j]]],
                               length(plcs),dim(reg_output)[2],byrow = TRUE)
  }
}

details = cbind(dataset[,1:(dim(dataset)[2]-1)],matrix(0,dim(dataset)[1],1),reg_output)
colnames(details) = c("L.neighbor","R.neighbor","codon",            
                      "codon.pos","gene","mat_peptide","stem_loop","CG","amino_acid",       
                      "output","ID","logLik","df","AIC","theta",
                      "NB.converged","rows_num","cols_num","P.GLR","P.test","P.logLik",
                      "P.df","P.AIC","P.converged","num_non_0_rows_1","output_sum") 

setwd("/Users/keren/Dropbox/covid/details/details_i_3_9/")
save(details, file = "missing_details_3_4")


##########3
# i=1
setwd("/Users/keren/Dropbox/covid/vars/datasets/reg_hashes/i_1/")
load("all_reg_hashes")
setwd("/Users/keren/Dropbox/covid/vars/datasets/syn_i_1/")
file_names = dir()
dataset = NULL
for (i in 1:length(file_names)){
  print(i)
  load(file_names[i])
  dataset = rbind(dataset,short_details[,1:11])
}
save(dataset,file = "dataset_i_1")

length_reg_hash = 15
reg_output = matrix(0,dim(dataset)[1],length_reg_hash)
reg_hash_names = names(all_reg_hash)
print(length(reg_hash_names))
for (j in 1:length(reg_hash_names)){
  print(j)
  plcs = which(dataset$reg_name2==reg_hash_names[j])
  if (length(plcs)>0){
    reg_output[plcs,] = matrix(all_reg_hash[[reg_hash_names[j]]],
                               length(plcs),dim(reg_output)[2],byrow = TRUE)
  }
}
details = cbind(dataset[,1:(dim(dataset)[2]-1)],matrix(0,dim(dataset)[1],1),reg_output)
colnames(details) = c("L.neighbor","R.neighbor","codon",            
                      "codon.pos","gene","mat_peptide","stem_loop","CG","amino_acid",       
                      "output","ID","logLik","df","AIC","theta",
                      "NB.converged","rows_num","cols_num","P.GLR","P.test","P.logLik",
                      "P.df","P.AIC","P.converged","num_non_0_rows_1","output_sum") 

setwd("/Users/keren/Dropbox/covid/details//")
save(details, file = "details_i_1")

details$output = 1
plcs = which(details$rows_num==0)
missing_dataset = dataset[plcs,]
setwd("/Users/keren/Dropbox/covid/vars/datasets/syn_i_1/")
save(missing_dataset,file = "missing_dataset_i_1")

setwd("/Users/keren/Dropbox/covid/vars/datasets/reg_hashes/i_1/")
load("missing_reg_hash_i_1")
setwd("/Users/keren/Dropbox/covid/vars/datasets/syn_i_1/")
load("missing_dataset_i_1")
dataset = missing_dataset
all_reg_hash = reg_hash
length_reg_hash = 15
reg_output = matrix(0,dim(dataset)[1],length_reg_hash)
reg_hash_names = names(all_reg_hash)
print(length(reg_hash_names))
for (j in 1:length(reg_hash_names)){
  print(j)
  plcs = which(dataset$reg_name2==reg_hash_names[j])
  if (length(plcs)>0){
    reg_output[plcs,] = matrix(all_reg_hash[[reg_hash_names[j]]],
                               length(plcs),dim(reg_output)[2],byrow = TRUE)
  }
}
details = cbind(dataset[,1:(dim(dataset)[2]-1)],matrix(0,dim(dataset)[1],1),reg_output)
colnames(details) = c("L.neighbor","R.neighbor","codon",            
                      "codon.pos","gene","mat_peptide","stem_loop","CG","amino_acid",       
                      "output","ID","logLik","df","AIC","theta",
                      "NB.converged","rows_num","cols_num","P.GLR","P.test","P.logLik",
                      "P.df","P.AIC","P.converged","num_non_0_rows_1","output_sum") 
details2 = details
setwd("/Users/keren/Dropbox/covid/details//")
load("details_i_1")
details$output = 1
plcs = which(details$rows_num==0)
details[plcs,] = details2
save(details, file = "details_i_1_all")


