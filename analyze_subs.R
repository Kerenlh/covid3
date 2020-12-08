#analyze subs:
setwd("/Users/keren/Dropbox/covid/vars/")
load("site_details")
load("all_subs")

subs_plcs = which(all_subs$before!="-" & all_subs$after!="-" & all_subs$uncertainty_flag==0)
subs = all_subs[subs_plcs,]
subs = data.frame(subs)
# subs_regions = matrix(0,dim(all_subs)[1],dim(regions)[2])
# subs_regions = data.frame(subs_regions)
# colnames(subs_regions) = colnames(regions)
subs_regions = NULL
for (i in 1:dim(subs)){
  print(i)
  sites = which(regions$site_pattern==subs$site_pattern[i])
  tmp = matrix(subs[i,],length(sites),dim(subs)[2],byrow = TRUE)
  subs_regions = rbind(subs_regions, cbind(tmp,regions[sites,]))
  # subs_regions [i,] = regions[which(regions$site_pattern==subs$site_pattern[i]),]
}
colnames(subs_regions) = c(colnames(subs),colnames(regions))
save(subs_regions,file = "subs_regions")

get_subs_mat = function(subs){
  subs = data.frame(subs)
  subs_mat = matrix(0,4,4)
  colnames(subs_mat) = c("A","C","G","T")
  rownames(subs_mat) = c("A","C","G","T")
  
  for (i in 1:dim(subs)[1]){
    print(i)
    subs_mat[subs$before[i][[1]],subs$after[i][[1]] ] = subs_mat[subs$before[i][[1]],subs$after[i][[1]] ]+1 
  }
  return(subs_mat)
}

all_subs_mat = get_subs_mat(subs_regions)
codon_subs_mat = get_subs_mat(subs_regions[subs_regions$gene!=0,])
mat_peptide_subs_mat = get_subs_mat(subs_regions[subs_regions$mat_peptide!=0,])
stem_loop_subs_mat = get_subs_mat(subs_regions[subs_regions$stem_loop!=0,])


