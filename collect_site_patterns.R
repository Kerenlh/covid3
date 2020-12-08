# setwd("/Users/keren/Dropbox/covid/vars")
# setwd("/Users/keren/Desktop/covid_files/missing_site_patterns/")
setwd("/Users/keren/Dropbox/covid/vars/missing_site_patterns2/")

file_names = dir()
# Missing files:
all_no_subs_plcs = NULL
all_site_patterns_plcs = NULL
all_site_patterns = NULL
for (i in 1:length(file_names)){
  tmp = strsplit(file_names[i],split="")[[1]]
  site_patterns = NULL
  no_subs_plcs = NULL
  site_patterns_plcs = NULL
  if (paste0(tmp[1:19],collapse="")=="site_patterns_plcs_"){
    site = as.numeric(paste0(tmp[20:length(tmp)],collapse = ""))
    load(file_names[i])
    load(paste0("site_patterns_",site,collaspe = ""))
    load(paste0("no_subs_plcs_",site,collapse = ""))
    all_site_patterns = c(all_site_patterns,site_patterns)
    all_site_patterns_plcs = c(all_site_patterns_plcs,site_patterns_plcs)
    all_no_subs_plcs = c(all_no_subs_plcs,no_subs_plcs)
   #  print(i)
   #  print(file_names[i])
   #  print(site_patterns_plcs)
   # # print(no_subs_plcs)
   #  print(site)
  }
}
all_no_subs_plcs = unique(all_no_subs_plcs)
all_site_patterns = unique(all_site_patterns)
all_site_patterns_plcs = unique(all_site_patterns_plcs)
# setwd("/Users/keren/Dropbox/covid/vars")
save(all_no_subs_plcs,file = "missing_all_no_subs_plcs")
save(all_site_patterns_plcs,file = "missing_all_site_patterns_plcs")
save(all_site_patterns,file = "missing_all_site_patterns")
setwd("/Users/keren/Dropbox/EC2_files/")
write.csv(all_site_patterns,file = "missing_all_site_patterns.csv")


file_names = dir()
all_no_subs_plcs = NULL
all_site_patterns_plcs = NULL
all_site_patterns = NULL
for (i in 1:length(file_names)){
  print(i)
  if (paste0(strsplit(file_names[i],split="")[[1]][1:12],collapse="")==
    "no_subs_plcs"){
    load(file_names[i])
    all_no_subs_plcs = c(all_no_subs_plcs,no_subs_plcs)
  }
  if (paste0(strsplit(file_names[i],split="")[[1]][1:14],collapse="")==
      "site_patterns_"){
    load(file_names[i])
    print(file_names[i])
    if (paste0(strsplit(file_names[i],split="")[[1]][15:19],collapse="")==
        "plcs_"){
      all_site_patterns_plcs = c(all_site_patterns_plcs,site_patterns_plcs)
    }else{
      all_site_patterns = c(all_site_patterns,site_patterns)
    }
  }
}


unique_site_patterns = unique(all_site_patterns)
all_letters = NULL
gaps = NULL
funny_letters = NULL
for (i in 1:length(unique_site_patterns)){
  # print(i)
  tmp = unique(strsplit(unique_site_patterns[i],split = "")[[1]])
  all_letters = c(all_letters, tmp)
  if (sum(is.element(tmp,"-"))>0){
    gaps = c(gaps,i)
  }
  if (sum(is.element(tmp,c("H", "N", "K", "Y", "B", "D", "R", "W" ,"V" ,"S" ,"M")))>0){
    funny_letters = c(funny_letters,i)
  }
}
saveRDS(all_site_patterns,file = "all_site_patterns.rds")
write.csv(all_site_patterns,file = "all_site_patterns.csv")
save(all_site_patterns_plcs, file = "all_site_patterns_plcs")
save(all_no_subs_plcs,file = "all_no_subs_plcs")
save(all_site_patterns,file = "all_site_patterns")