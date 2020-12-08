setwd("/Users/keren/Dropbox/covid/vars/datasets/base/")

file_names = dir()
dataset = NULL
count = 0
for (i in 1:length(file_names)){
  print(i)
  if (paste0(strsplit(file_names[i],split = "")[[1]][1:14]
             ,collapse = "")=="short_details."){
    load(file_names[i])
    count = count + dim(short_details)[1]
    plcs = which(duplicated(short_details$reg_name2)==FALSE)
    dataset = rbind(dataset,short_details[plcs,])
  }
}

plcs = which(duplicated(dataset$reg_name2)==FALSE)
dataset = dataset[plcs,]
setwd("/Users/keren/Dropbox/EC2_files/datasets/")
save(dataset,file = "unique_dataset_base")