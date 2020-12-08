# Apply calc_missing_likelihood_AIC_NA before!

# setwd("/Users/keren/Dropbox/covid/details/details_i_3_9/")
# load("missing_details_6_9")
# details2 = details
# load("details_6_9_missing_likelihood")
# details = rbind(details,details2)
# rm(details2)

setwd("/Users/keren/Dropbox/covid/details/")
load("details_non_syn")
details2 = details
load("details_i_1_all_missing_likelihood")
details = rbind(details,details2)
rm(details2)

setwd("/Users/keren/Dropbox/covid/vars/")
load("iterate_vals_codons")
load("codons_table3")
data = codons_table
source("/Users/keren/Dropbox/covid/debug_functions.R")

##############
# Functions: #
##############
get_model_ids = function(details,iterate_vals){
  colnames_model_ids = c(colnames(iterate_vals),"output","ID")
  model_ids = matrix(1,dim(details)[1],length(colnames_model_ids))
  colnames(model_ids) = colnames_model_ids
  
  model_ids[which(details$codon.pos==4),"codon.pos"] = 2
  model_ids[which(is.na(details$codon.pos)),"codon.pos"] = 3
  model_ids[,c("R.neighbor","L.neighbor")] = 0
  # 1 = divide according to neighbors with respect to codon.pos
  model_ids[which(details$codon.pos==1 & details$L.neighbor<5 & is.na(details$R.neighbor)),c("R.neighbor","L.neighbor")] = 1 
  model_ids[which(details$codon.pos==3 & is.na(details$L.neighbor) & details$R.neighbor<5),c("R.neighbor","L.neighbor")] = 1
  # 2 = have neighbors as an explaining variable outside the codon 
  model_ids[which(details$codon.pos==1 & details$L.neighbor==5 & is.na(details$R.neighbor)),c("R.neighbor","L.neighbor")] = 2 
  model_ids[which(details$codon.pos==3 & is.na(details$L.neighbor) & details$R.neighbor==5),c("R.neighbor","L.neighbor")] = 2
  # 3 = don't include neighbors
  model_ids[which(details$codon.pos==1 & is.na(details$L.neighbor) & is.na(details$R.neighbor)),c("R.neighbor","L.neighbor")] = 3 
  model_ids[which(details$codon.pos==3 & is.na(details$L.neighbor) & is.na(details$R.neighbor)),c("R.neighbor","L.neighbor")] = 3
  # codon.pos==2 is the same for all options and remaining codon.pos are foreign for different options
  model_ids[which(details$codon.pos==2 & is.na(details$L.neighbor) & is.na(details$R.neighbor)),c("R.neighbor","L.neighbor")] = 5
  
  # codon.pos==4:
  model_ids[which((details$codon.pos==4 | is.na(details$codon.pos)) & details$L.neighbor<5),c("L.neighbor")] = 1
  model_ids[which((details$codon.pos==4 | is.na(details$codon.pos)) & details$L.neighbor==5),c("L.neighbor")] = 2
  model_ids[which((details$codon.pos==4 | is.na(details$codon.pos)) & is.na(details$L.neighbor)),c("L.neighbor")] = 3
  model_ids[which((details$codon.pos==4 | is.na(details$codon.pos)) & details$R.neighbor<5),c("R.neighbor")] = 1
  model_ids[which((details$codon.pos==4 | is.na(details$codon.pos)) & details$R.neighbor==5),c("R.neighbor")] = 2
  model_ids[which((details$codon.pos==4 | is.na(details$codon.pos)) & is.na(details$R.neighbor)),c("R.neighbor")] = 3
  
  
  model_ids[which(details$mat_peptide==3),"mat_peptide"] = 2
  model_ids[which(is.na(details$mat_peptide)),"mat_peptide"] = 3
  model_ids[which(details$CG==4),"CG"] = 2
  model_ids[which(is.na(details$CG)),"CG"] = 3
  model_ids[which(details$codon==65),"codon"] = 2
  model_ids[which(is.na(details$codon)),"codon"] = 3
  model_ids[which(details$amino_acid==22),"amino_acid"] = 2
  model_ids[which(is.na(details$amino_acid)),"amino_acid"] = 3
  model_ids[which(details$gene==12),"gene"] = 2
  model_ids[which(is.na(details$gene)),"gene"] = 3
  model_ids[which(details$stem_loop==4),"stem_loop"] = 2
  model_ids[which(is.na(details$stem_loop)),"stem_loop"] = 3
  model_ids[which(details$output==1 | details$output==2),"output"] = 1
  model_ids[which(details$output==3 | details$output==4),"output"] = 2
  model_ids[which(details$output==5),"output"] = 3
  model_ids[which(details$output==6 | details$output==7 | 
                    details$output==8 | details$output==9),"output"] = 4
  
  
  # model_ids[which(details$base_backwords==5),"base_backwords"] = 2
  # model_ids[is.na(details$base_backwords),"base_backwords"] = 3
  
  model_ids[,"ID"] = as.numeric(apply(model_ids[,1:(dim(model_ids)[2]-1)],1,paste0,collapse = ""))
  model_ids = data.frame(model_ids)
  details_ids = apply(details[,1:(dim(model_ids)[2]-1)],1,paste0,collapse = "")
  tmp = list(model_ids,details_ids)
  return(tmp)
}

find_id_plcs = function(model_ids,unique_models,unique_model_nums_output,unique_model_nums_output2,model_details_ids,
                        model_nums_details_ids2,model_nums2){
  for (i in 1:dim(unique_models)[1]){
    print(i)
    curr_id = unique_models$ID[i]
    plcs = which(model_ids$ID==curr_id)
    output_pos = dim(model_ids)[2]-1
    if(substring(curr_id,first = 1,last = 1)=="1" & substring(curr_id,first = 6,last = 7)=="55"){next}
    if(substring(curr_id,first = 1,last = 1)=="1" & substring(curr_id,first = 6,last = 7)!="55"){
      curr_id2 = as.numeric(paste0(substring(curr_id,first = 1,last=5),"55",substring(curr_id,first = 8,last=output_pos),collapse = ""))
      plcs = c(plcs, which(model_ids$ID==curr_id2))
    }
    unique_models[i,colnames(model_ids)] = model_ids[plcs[1],]
    unique_models$num_models[i] = length(plcs)
    unique_models$sum_rows[i] = sum(details[plcs,"rows_num"])
    unique_models$logLik.NB[i] = sum(details[plcs,"logLik"])
    unique_models$logLik.P[i] = sum(details[plcs,"P.logLik"])
    
    if (substring(curr_id,first = output_pos,last = output_pos)==3){
      unique_models$num_missing_rows[i] = dim(data)[1]-unique_models$sum_rows[i]
    }
    if (substring(curr_id,first = output_pos,last = output_pos)==2){
      unique_models$num_missing_rows[i] = 2*dim(data)[1]-unique_models$sum_rows[i]
    }
    syn_non_syn_rows_num = length(which(data$syn_exposure>0))+length(which(data$non_syn_exposure>0))
    if (substring(curr_id,first = output_pos,last = output_pos)==1){ # syn/non_syn together have 27144 rows
      unique_models$num_missing_rows[i] = syn_non_syn_rows_num-unique_models$sum_rows[i]
    }
    if (substring(curr_id,first = output_pos,last = output_pos)==4){ # ACGT
      unique_models$num_missing_rows[i] = 4*dim(data)[1]-unique_models$sum_rows[i]
    }
    
    unique_models$total_num_models[i] = unique_models$num_models[i]+unique_models$num_missing_rows[i]
    
    unique_models$num_missing_models[i] =unique_models$num_missing_rows[i] # assumes all missing models include only 1 row!
    
    unique_models$df.NB[i] = sum(details$df[plcs])+unique_models$num_missing_models[i]
    unique_models$df.P[i] = sum(details$P.df[plcs])+unique_models$num_missing_models[i]
    
    unique_models$AIC.NB[i] = 2*(unique_models$df.NB[i]-unique_models$logLik.NB[i])
    unique_models$AIC.P[i] = 2*(unique_models$df.P[i]-unique_models$logLik.P[i])
    unique_models$NA.NB[i] = sum(details$added_log_lik[plcs])
    unique_models$NA.P[i] = sum(details$added_log_lik.P[plcs])
    
    if (length(which(details$added_log_lik[plcs]==1))>0){
      good_plcs.NB = plcs[-which(details$added_log_lik[plcs]==1)]
    }else{
      good_plcs.NB = plcs
    }
    unique_models[i,26:31] = summary(details[good_plcs.NB,"theta"])
    unique_models[i,32] = weighted.mean(details[good_plcs.NB,"theta"],
                                        details[good_plcs.NB,"rows_num"]/sum(details[good_plcs.NB,"rows_num"]))
    unique_models[i,33] = sum(details$output_0_or_1_row[plcs])
  }
  return(unique_models)
}



########

tmp = get_model_ids(details,iterate_vals)
model_ids = tmp[[1]]; model_details_ids = tmp[[2]]
tmp = unique(model_ids[,"ID"])

colnames_unique_models = c(colnames(model_ids),"num_models","logLik.NB","df.NB",
                           "num_missing_models","AIC.NB","NA.NB","logLik.P","df.P",
                           "AIC.P","NA.P","sum_rows","total_num_models","num_missing_rows",
                           "theta_min","theta_1st Qu.","theta_Median","theta_Mean","theta3rd Qu.",
                           "theta_Max","theta_weighted_mean","output_0_or_1_row_models")
unique_models = matrix(0,length(tmp),length(colnames_unique_models))
colnames(unique_models) = colnames_unique_models
unique_models = data.frame(unique_models)
unique_models$ID = tmp

unique_models = find_id_plcs(model_ids,unique_models,unique_model_nums_output,unique_model_nums_output2,model_details_ids,
                             model_nums_details_ids2,model_nums2)
unique_models = unique_models[-which(unique_models$codon.pos==0),]
setwd("/Users/keren/Dropbox/covid/models")
save(unique_models,file = "unique_models_i_1_2")

# load("unique_models_base")
# load("unique_models_base_no_repetitions") # corrected in compare_results to avoid repetitions.
plcs = which(unique_models$num_missing_models!=0)
unique_models2 = unique_models[plcs,]
plcs = which(unique_models$num_missing_rows==0)
tmp2 = unique_models[order(unique_models$AIC.NB,decreasing = FALSE),]

plot(tmp2$AIC.NB, main = "AIC score for different models",ylab = "AIC")
points(tmp2$AIC.P
       ,col = "blue")
legend(300,90000, c("NB AIC", "Poisson AIC"), lty=c(1,1), lwd=c(2.5,2.5),col=c("black","blue")) 
write.csv(tmp2,file = "unique_models_base.csv")

tmp3 = unique_models[order(unique_models$AIC.P),]
plot(tmp3$AIC.P)

tmp4 = unique_models[order(apply(unique_models[,c("AIC.NB","AIC.P")],1,min)),]
write.csv(tmp4,file = "unique_models_base_tmp4.csv")
plot(apply(tmp4[,c("AIC.NB","AIC.P")],1,min))
points(tmp4$AIC.P,col = "red")
points(tmp4$AIC.NB,col = "blue")

plcs = which(tmp4$AIC.P-tmp4$AIC.NB<500)

# Debug:

# i =3964668 # i is the relevant row number in details (i's:184019,371386 )
i = 3310023
tmp = get_data(i,iterate_vals)
curr_data = tmp[[1]]; curr_details = tmp[[2]]
# tmp = prepare_data_rm_cols(curr_data)
# curr_data = tmp[[1]]; curr_model_names = tmp[[2]]
# model.nb = find_model.nb.data(curr_data,curr_model_names)
# model.p = find_model.P.data(curr_data, curr_model_names)


# tmp = which(is.na(details$codon.pos) & details$L.neighbor<5 & details$R.neighbor<5 & is.na(details$codon) &
#               details$gene<14 & is.na(details$stem_loop) & 
#               is.na(details$mat_peptide) & is.na(details$CG) & details$output==5 & 
#               details$amino_acid==22 &
#               details$base_backwords==5)
# load("unique_models_base_missing_models")
# missing_models = unique_models
# load("unique_models_base")
# tmp = is.element(tmp2,unique_models$ID)
# 
# tmp2 = get_model_ids(missing_models,iterate_vals)
# model_ids2 = tmp2[[1]]; model_details_ids2 = tmp2[[2]]
# tmp2 = unique(model_ids2[,"ID"])
# save(redundant_ids,file = "redundant_ids")
# tmp = which(missing_models)
# 
# tmp = is.element(model_ids2$ID,redundant_ids)
# missing_models_non_redundant = missing_models[which(tmp==FALSE),]
# tmp3 = get_model_ids(missing_models_non_redundant,iterate_vals)
# model_ids3 = tmp3[[1]]; model_details_ids3 = tmp2[[2]]
# tmp3 = unique(model_ids3[,"ID"])
# save(missing_models_non_redundant,file = "missing_models_non_redundant")

