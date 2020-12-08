setwd("/Users/keren/Dropbox/covid/details/All_details_non_syn/")

all_details = NULL
nums = NULL

file_names = dir()
for (i in c(1:length(file_names))){
  if(paste0(strsplit(file_names[i],split = "")[[1]][1:8],collapse = "")=="details."){
    # print(file_names[i])
    print(i)
    load(file_names[i])
    tmp = strsplit(file_names[i],split = "\\.")[[1]]
    nums = c(nums,tmp[length(tmp)])
    if (is.null(dim(details))){
      print(paste(file_names[i],"is NULL"))
    }else{
      # print(details[1:10,])
      all_details = rbind(all_details,details)
    }
  }
}

nums = as.numeric(nums)
missing_nums = which(is.element(1:max(nums),nums)==FALSE)
  
details = data.frame(all_details)
tmp = which(is.na(details$theta))
print(paste("theta=NA: ",length(tmp)))
print(table(details$output[tmp]))
# details = details[-tmp,]
tmp = which(is.na(details$P.GLR))
print(paste("GLR.nb.P=NA: ",length(tmp),table(details$output[tmp])))
# details = details[-tmp,]
tmp = which(details$rows_num==1)
print(paste("rows_num=0: ",length(tmp),table(details$output[tmp])))
tmp = which(is.na(details$logLik))
print(paste("logLik=NA: ",length(tmp),table(details$output[tmp])))
# print(dim(details))
tmp = which(details$NB.converged==0)
print(paste("NB.converged=0: ",length(tmp),table(details$output[tmp])))
# print(dim(details))
# details = details[-tmp,]
tmp = which(details$P.converged==0)
print(paste("P.converged=0: ",length(tmp),table(details$output[tmp])))
# print(dim(details))


setwd("/Users/keren/Dropbox/covid/details/")
save(details,file = "details_non_syn")
# setwd("C:\\Users\\user\\Dropbox\\mtDNA_tree_Build_17")
# load("details_rRNA")
# details = rbind(details,details2)
# save(details,file = "details_control_base2")


# # setwd("C:\\Users\\Keren\\Dropbox\\mtDNA_tree_Build_17\\num_models")
# setwd("C:\\Users\\user\\Dropbox\\mtDNA_tree_Build_17\\num_models")
# file_names = dir()
# all_models = NULL
# for (i in c(1:length(file_names))){
#   print(i)
#   load(file_names[i])
#   if (is.null(dim(details))){
#     print(paste(file_names[i],"is NULL"))
#   }else{
#     all_models = rbind(all_models,details)
#   }
# }
# 
# model_nums = data.frame(all_models)
# # model_nums = model_nums[-which(model_nums$rows_num==0),]
# for (i in 1:13){
#   model_nums[,i] = as.numeric(as.character(model_nums[,i])) 
# }
# model_nums2 = model_nums
# setwd("C:\\Users\\user\\Dropbox\\mtDNA_tree_Build_17")
# save(model_nums2,file = "model_nums2")

