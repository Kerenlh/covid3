setwd("/Users/keren/Dropbox/covid/details/")
# load("details_3_4")
load("details_i_1_all")
plcs = which(details$cols_num==0)
if(length(plcs)>0){
  details = details[-plcs,]
}
details$output = 1

setwd("/Users/keren/Dropbox/covid/vars//")
load("iterate_vals_codons")
load("codons_table3")

source("/Users/keren/Dropbox/covid/debug_functions.R")

added_log_lik = matrix(0,dim(details)[1],1)
added_log_lik.P = matrix(0,dim(details)[1],1)
output_0_or_1_row = matrix(0,dim(details)[1],1)
details = cbind(details,added_log_lik,added_log_lik.P, output_0_or_1_row)

plcs = which(is.na(details$AIC))
if (length(plcs)>0){
  log_liks = matrix(0,length(plcs),1)
  for (i in 1:length(plcs)){
    print(i)
    tmp = get_data(details[plcs[i],],iterate_vals)
    curr_data = tmp[[1]];# curr_details = tmp[[2]]
    ys = curr_data[,dim(curr_data)[2]]
    curr_lik = 1
    for (j in 1:length(ys)){
      ys_j = ys[j]
      # curr_lik = curr_lik*exp(-ys_j)*(ys_j^ys_j)/factorial(ys_j)
      curr_lik = curr_lik*dpois(ys_j,ys_j)
      print(curr_lik)
      print(j)
      if(curr_lik==0){
        break
      }
    }
    log_liks[i,1] = log(curr_lik)
  } 
  details[plcs,"logLik"] = log_liks
  # added_log_lik = matrix(0,dim(details)[1],1)
  added_log_lik[plcs] = 1
  # details = cbind(details,added_log_lik)
  # details$df[plcs] = details$num_non_0_rows_1[plcs]
  details$df[plcs] = details$rows_num[plcs]
  # save(details,file = "details_non_codons_missing_models")
}


#Poisson correction:
plcs = which(is.na(details$P.AIC) | (details$P.converged==0 & 
                                       (details$P.df==0 | (details$P.df!=0 & details$P.AIC>1e5)) ))
if (length(plcs)>0){
  log_liks = matrix(0,length(plcs),1)
  for (i in 1:length(plcs)){
    print(i)
    tmp = get_data(details[plcs[i],],iterate_vals)
    curr_data = tmp[[1]];# curr_details = tmp[[2]]
    ys = curr_data[,dim(curr_data)[2]]
    curr_lik = 1
    for (j in 1:length(ys)){
      ys_j = ys[j]
      # curr_lik = curr_lik*exp(-ys_j)*(ys_j^ys_j)/factorial(ys_j)
      curr_lik = curr_lik*dpois(ys_j,ys_j)
      if(curr_lik==0){
        break
      }
    }
    log_liks[i,1] = log(curr_lik)
  }
  details[plcs,"P.logLik"] = log_liks
  
  # added_log_lik.P = matrix(0,dim(details)[1],1)
  added_log_lik.P[plcs] = 1
  # details = cbind(details,added_log_lik.P)
  
  # details$P.df[plcs] = details$num_non_0_rows_1[plcs]
  details$P.df[plcs] = details$rows_num[plcs]
  # save(details,file = "details_non_codons_missing_models")
}

#models with 1 row and output!=0
plcs  = which(is.na(details$NB.converged) & is.na(details$AIC)==FALSE 
              & details$output_sum!=0)
if (length(plcs)>0){
  details$P.df[plcs] = 1
  details$df[plcs] = 1
  vals = unique(details$output_sum[plcs])
  for (i in 1:length(vals)){
    curr_val = vals[i]
    curr_plcs = which(is.na(details$NB.converged) & is.na(details$AIC)==FALSE
                      & details$output_sum==curr_val)
    curr_lik = exp(-curr_val)*(curr_val^curr_val)/factorial(curr_val)
    print(log(curr_lik))
    print(length(curr_plcs))
    details[curr_plcs,c("logLik","P.logLik")] = log(curr_lik)
    details[curr_plcs,c("added_log_lik","added_log_lik.P")] = 1
  }
  
  # save(details,file = "details_non_codons_missing_models")
  
  plcs = which(is.na(details$NB.converged) & is.na(details$AIC)==FALSE)
  details$added_log_lik[plcs] = 1
  details$added_log_lik.P[plcs] = 1
  # save(details,file = "details_non_codons_missing_models")
}

plcs = which(details$rows_num==1 | details$output_sum==0)
if (length(plcs)>0){
  # output_0_or_1_row = matrix(0,dim(details)[1],1)
  output_0_or_1_row[plcs] = 1
  # details = cbind(details,output_0_or_1_row)
}
setwd("/Users/keren/Dropbox/covid/details/")
save(details,file = "details_i_1_all_missing_likelihood")
