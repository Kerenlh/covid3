setwd("/a/home/cc/math/kerenlh/EC2_files/Runs/unique_regs/base")
wait = runif(1, 0, 30)
Sys.sleep(wait)
if (file.exists("begin")){
  load("begin")
}else{
  begin = 1
}

old_begin = begin

job_length = 500#10000 #last-begin+1
last = begin+job_length-1
begin = last+1
save(begin,file = "begin")

begin = old_begin
print(begin)
print(last)

setwd("/a/home/cc/math/kerenlh/EC2_files")
# setwd("/Users/keren/Dropbox/EC2_files/datasets/")
# load("unique_dataset_i_3_9")
load("unique_dataset_base")
load("codons_table3")
load("iterate_vals")
# source("/Users/keren/Dropbox/covid/debug_functions.R")
source("/a/home/cc/math/kerenlh/EC2_files/debug_functions.R")
library(hash)
library(digest)

reg_hash = hash()

for (i in begin:last){
  print(i)
  tmp = get_data(dataset[i,],iterate_vals = iterate_vals)
  curr_data = tmp[[1]];
  curr_details = tmp[[2]]
  ys = curr_data[,dim(curr_data)[2]]
  added_log_lik = added_log_lik.P = output_0_or_1_row = 0
  curr_details = cbind(curr_details,added_log_lik,added_log_lik.P, output_0_or_1_row)
  if (is.na(curr_details[,"AIC"])){
    curr_lik = 1
    for (j in 1:length(ys)){
      ys_j = ys[j]
      # curr_lik = curr_lik*exp(-ys_j)*(ys_j^ys_j)/factorial(ys_j)
      curr_lik = curr_lik*dpois(ys_j,ys_j)
      if(curr_lik==0){
        break
      }
    }
    curr_details[,"logLik"] = log(curr_lik)
    # added_log_lik = matrix(0,dim(details)[1],1)
    added_log_lik = 1
    curr_details[,"df"] = curr_details[,"rows_num"]
  }
  
  #Poisson correction:
  if (is.na(curr_details[,"P.AIC"]) | (curr_details[,"P.converged"]==0 &
                                       (curr_details[,"P.df"]==0 | (curr_details[,"P.df"]!=0 & curr_details[,"P.AIC"]>1e5)))){
    curr_lik = 1
    for (j in 1:length(ys)){
      ys_j = ys[j]
      curr_lik = curr_lik*dpois(ys_j,ys_j)
      if(curr_lik==0){
        break
      }
    }
    curr_details[,"P.logLik"] = log(curr_lik)
    added_log_lik.P = 1
    curr_details[,"P.df"] = curr_details[,"rows_num"] 
  }
  
  #models with 1 row and output!=0
  if (is.na(curr_details[,"NB.converged"]) & is.na(curr_details[,"AIC"])==FALSE 
      & curr_details[,"output_sum"]!=0){
    curr_details[,"P.df"] = 1
    curr_details[,"df"] = 1
    curr_val = curr_details[,"output_sum"]
    curr_lik = exp(-curr_val)*(curr_val^curr_val)/factorial(curr_val)
    if (is.na(curr_details[,"NB.converged"]) & is.na(curr_details[,"AIC"])==FALSE
        & curr_details[,"output_sum"]==curr_val){
      curr_details[,c("logLik","P.logLik")] = log(curr_lik)
      curr_details[,c("added_log_lik","added_log_lik.P")] = 1
    }
    
  }
  if(is.na(curr_details[,"NB.converged"]) & is.na(curr_details[,"AIC"])==FALSE){
    curr_details[,"added_log_lik"] = 1
    curr_details[,"added_log_lik.P"] = 1
  }
  
  if(curr_details[,"rows_num"]==1 | curr_details[,"output_sum"]==0){
    output_0_or_1_row = 1
  }
  reg_hash[[dataset$reg_name2[i]]] = 
    curr_details[(length(iterate_vals)+3):length(curr_details)]
}

setwd("/a/home/cc/math/kerenlh/EC2_files/Runs/unique_regs/base")
save(reg_hash,file = paste0("reg_hash.",begin,collapse = ""))
stop("so that the workspace is not saved")

