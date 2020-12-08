setwd("/Users/keren/Desktop/covid_files//trees")
file_names = dir()
tmp = which(file_names=="subs")
file_names = file_names[-tmp]
site_patterns = NULL
for (i in 1:length(file_names)){
  #print(i)
  tmp = strsplit(file_names[i],split = "")[[1]]
  site_pattern = as.numeric(paste0(tmp[11:length(tmp)],collapse=""))
  site_patterns = c(site_patterns,site_pattern)
}
tmp = which(is.element(c(1:13350),site_patterns)==FALSE)
missing_plcs = tmp


info = file.info(file_names)
empty_files_plcs = which(info$size<500000)
empty_files = NULL
for (i in 1:length(empty_files_plcs)){
  print(i)
  tmp = strsplit(file_names[empty_files_plcs[i]],split = "")[[1]]
  empty_file = as.numeric(paste0(tmp[11:length(tmp)],collapse=""))
  empty_files = c(empty_files,empty_file)
}


setwd("/Users/keren/Dropbox/EC2_files/")
write.csv(missing_plcs,file = "missing_plcs.csv")
write.csv(empty_files, file = "empty_files.csv")

load("subs")