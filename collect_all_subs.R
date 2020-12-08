setwd("/Users/keren/Desktop/covid_files/subs/")
file_names = dir()

# all_subs = NULL
# nums = NULL
# for (i in 1:length(file_names)){
#   load(file_names[i])
#   print(i)
#   tmp = strsplit(file_names[i],split = "")[[1]]
#   nums = c(nums,as.numeric(paste0(tmp[6:length(tmp)],collapse = "")))
#   all_subs = rbind(all_subs,subs)
# }

all_subs = data.frame(all_subs)
# table(as.numeric(all_subs$site_pattern))
missing_subs = which(is.element(1:13366,nums)==FALSE)
setwd("/Users/keren/Dropbox/covid/vars/")
save(all_subs,file = "all_subs")

