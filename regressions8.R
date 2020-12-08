create_details = function(parallel_num){
  # setwd("/Users/keren/Dropbox/mtDNA_tree_Build_17")
  # data.codon_nbs = read.csv("codon_exposure_neighbors_data_nbs.csv",string=F
  
  # setwd("/Users/keren/Dropbox/covid/vars/")
  setwd("/home/ubuntu/Dropbox/EC2_files/")
  load("codons_table2")
  plcs = which(is.na(codons_table$non_syn))
  codons_table = codons_table[-plcs,]
  plcs = which(codons_table$exposure==0)
  codons_table = codons_table[-plcs,]
  
  
  library(MASS)
  ##############
  # Contrasts: #
  ##############
  create_contrasts = function(k){
    c<-contr.treatment(k)
    my.coding<-matrix(rep(1/k, k*(k-1)), ncol=(k-1))
    my.simple<-c-my.coding
    return(my.simple)
  }
  
  add_contrasts = function(data,cols){
    for (i in 1:length(cols)){
      print(cols[i])
      data[,cols[i]] = as.factor(data[,cols[i]])
      contrasts(data[,cols[i]]) = create_contrasts(length(table(data[,cols[i]])))
    }
    return(data)
  }
  codons_table$mat_peptide[which(codons_table$mat_peptide!="0")] = 1
  codons_table$mat_peptide = as.numeric(codons_table$mat_peptide)
  tmp = c("L.neighbor","R.neighbor","codon","codon.pos","gene",            
          "stem_loop","amino_acid","base")
  codons_table = add_contrasts(codons_table,tmp)
  # tmp = c("codon.pos","regions","regions2","base","before","codon.a","codon.b","codon.c",
  #         "codon","codon_neighbor","before","codon.a_backwords","codon.b_backwords",
  #         "codon.c_backwords", "R.neighbor_backwords","L.neighbor_backwords","R.neighbor",
  #         "L.neighbor","codon","amino_acid","amino_acid_backwords","CG_position")
  # data.codon_nbs[data.codon_nbs$regions2==0,"regions2"]="none"
  # data.codon_nbs = add_contrasts(data.codon_nbs,tmp)
  
  ##############
  # functions: #
  ##############
  
  nb.glm.theta.data = function(theta,data,curr_model_names){
    # print(theta)
    model.nb = glm(curr_model_names,
                   family = negative.binomial(theta),data = data,maxit = 1000)
    logLik_model.nb = logLik(model.nb)
    return(logLik_model.nb)
  }
  
  find_model.nb.data = function(curr_data,curr_model_names){
    # if (colnames(curr_data)[dim(curr_data)[2]]=="syn"){
    #   curr_model_names = as.formula(paste(colnames(curr_data)[dim(curr_data)[2]],"~",
    #                                       ".-exposure-syn_exposure+offset(log(exposure)+log(syn_exposure))"))
    #   curr_data = curr_data[,-which(colnames(curr_data)=="non_syn_exposure")]
    # }else if (colnames(curr_data)[dim(curr_data)[2]]=="non_syn"){
    #   curr_model_names = as.formula(paste(colnames(curr_data)[dim(curr_data)[2]],"~",
    #                                       ".-exposure-non_syn_exposure+offset(log(exposure)+log(non_syn_exposure))"))
    #   curr_data = curr_data[,-which(colnames(curr_data)=="syn_exposure")]
    # }else{
    #   curr_data = curr_data[,-c(which(colnames(curr_data)=="syn_exposure"),
    #                             which(colnames(curr_data)=="non_syn_exposure"))]
    #   curr_model_names = as.formula(paste(colnames(curr_data)[dim(curr_data)[2]],"~",
    #                                       " . -exposure + offset(log(exposure))"))
    # }
    
    theta <- tryCatch({
      a = optimize(nb.glm.theta.data, data=curr_data,curr_model_names = curr_model_names,
                   interval=c(0.05, 100), maximum=TRUE)
      theta = as.numeric(a[1])
    }, warning = function(war) {
      a = optimize(nb.glm.theta.data, data=curr_data,curr_model_names = curr_model_names,
                   interval=c(0.05, 100), maximum=TRUE)
      theta = as.numeric(a[1])
      return(theta)
    }, error = function(err) {
      print(paste("MY_ERROR find_model.nb.data:  ",err))
      theta = -1
      return(theta)
    }, finally = {
    }) # END tryCatch
    
    # a = optimize(nb.glm.theta.data, data=curr_data, interval=c(0.05, 100), maximum=TRUE)
    # theta = as.numeric(a[1])
    if(theta!=-1){
      model.nb = glm(curr_model_names,
                     family = negative.binomial(theta),data = curr_data,maxit = 1000)
    }else{model.nb = 0}
    
    return(model.nb)
  }
  
  find_model.P.data = function(data,curr_model_names){
    # if (colnames(data)[dim(data)[2]]=="syn"){
    #   data = data[,-which(colnames(data)=="non_syn_exposure")] 
    #   curr_model_names = as.formula(paste(colnames(data)[dim(data)[2]],"~",
    #                                       ".-exposure-syn_exposure+offset(log(exposure)+log(syn_exposure))"))
    # }else if(colnames(data)[dim(data)[2]]=="non_syn") {
    #   data = data[,-which(colnames(data)=="syn_exposure")]
    #   curr_model_names = as.formula(paste(colnames(data)[dim(data)[2]],"~",
    #                                       ".-exposure-non_syn_exposure+offset(log(exposure)+log(non_syn_exposure))"))
    # }
    # else{
    #   data = data[,-c(which(colnames(data)=="syn_exposure"),which(colnames(data)=="non_syn_exposure"))]
    #   curr_model_names = as.formula(paste(colnames(data)[dim(data)[2]],"~",
    #                                       ".-exposure+offset(log(exposure))"))
    # }
    
    model.p <- tryCatch({
      model.p = glm(curr_model_names,data = data,family = "poisson")
      # GLR = 2*(logLik(model.nb)[1]-logLik(model.p)[1])
    }, warning = function(war) {
      model.p = glm(curr_model_names,data = data,family = "poisson")
      # GLR = 2*(logLik(model.nb)[1]-logLik(model.p)[1])
      return(model.p)
    }, error = function(err) {
      print(paste("MY_ERROR find_model.P.data:  ",err))
      model.p = -1
      return(model.p)
    }, finally = {
    }) # END tryCatch
    
    return(model.p)
  }
  
  remove_linear_dependant_cols = function(curr_data,log_plcs){
    #recursive!
    rankifremoved <- sapply(1:ncol(curr_data), function (x) qr(curr_data[,-x])$rank)
    if (length(unique(rankifremoved))>1){
      col_to_remove = which(rankifremoved == max(rankifremoved))
      col_to_remove = col_to_remove[1]   #removes only 1 col each time and makes sure it's not the output
      if(col_to_remove>=(ncol(curr_data)-length(log_plcs)) ){
        return(curr_data) # doesn't remove output and exposure cols
      }else{
        curr_data = curr_data[,-col_to_remove]
        remove_linear_dependant_cols(curr_data,log_plcs)
      }
    }else{
      return(curr_data) 
    }
  }
  
  prepare_data_rm_cols = function(curr_data){
    if (colnames(curr_data)[dim(curr_data)[2]]=="syn" & length(unique(curr_data$syn_exposure))>1 ){
      curr_model_names = as.formula(paste(colnames(curr_data)[dim(curr_data)[2]],"~",
                                          ".-exposure-syn_exposure+offset(log(exposure)+log(syn_exposure))"))
      # curr_data = curr_data[,-which(colnames(curr_data)=="non_syn_exposure")]
      log_plcs = c(which(colnames(curr_data)=="syn_exposure"), which(colnames(curr_data)=="exposure"))
    }else if (colnames(curr_data)[dim(curr_data)[2]]=="non_syn" & length(unique(curr_data$non_syn_exposure))>1){
      curr_model_names = as.formula(paste(colnames(curr_data)[dim(curr_data)[2]],"~",
                                          ".-exposure-non_syn_exposure+offset(log(exposure)+log(non_syn_exposure))"))
      # curr_data = curr_data[,-which(colnames(curr_data)=="syn_exposure")]
      log_plcs = c(which(colnames(curr_data)=="non_syn_exposure"), which(colnames(curr_data)=="exposure"))
    }else{
      # curr_data = curr_data[,-c(which(colnames(curr_data)=="syn_exposure"),
      # which(colnames(curr_data)=="non_syn_exposure"))]
      curr_model_names = as.formula(paste(colnames(curr_data)[dim(curr_data)[2]],"~",
                                          " . -exposure + offset(log(exposure))"))
      log_plcs = which(colnames(curr_data)=="exposure")
    }
    if(dim(curr_data)[2]>4){
      curr_data_input = curr_data[,1:(ncol(curr_data)-4),drop=FALSE ] # remove exposure, syn_exposure, non_syn_exposure and output
      factor_cols = which(sapply(curr_data_input, class)=="factor")
      for (i in factor_cols){
        curr_data_input[,i] = as.numeric(factor(curr_data_input[,i]))
      }
      curr_data_factored = cbind(curr_data_input,log(curr_data[,log_plcs]),curr_data[,ncol(curr_data)])
    }else{
      curr_data_factored = cbind(log(curr_data[,log_plcs]),curr_data[,ncol(curr_data)])
    }
    
    colnames(curr_data_factored)[(ncol(curr_data)-3):ncol(curr_data_factored)] = colnames(curr_data)[c(log_plcs,ncol(curr_data))]
    # if((ncol(curr_data)-4)==1){
    #   colnames(curr_data_factored)[1] = colnames(curr_data)[1]
    # }
    curr_data_factored = remove_linear_dependant_cols(curr_data_factored,log_plcs)
    curr_data = curr_data[,colnames(curr_data_factored)]
    
    tmp = list(curr_data,curr_model_names)
    return(tmp)
  }
  
  find_df = function(model){
    df = as.numeric(0.5*model$aic+logLik(model))
    return(df)
  }
  
  need_to_divide.GLR = function(model_nums,details,num_missing_models){
    sep_logLik = 0
    sep_df = num_missing_models
    model_nums = as.numeric(model_nums)
    for (r in model_nums[1:(length(model_nums)-1)]){
      sep_logLik = sep_logLik + as.numeric(details[r,"logLik"])
      sep_df = sep_df + as.numeric(details[r,"df"])
    }
    GLR = 2*(sep_logLik-as.numeric(details[model_nums[length(model_nums)],"logLik"]))
    df_diff = sep_df-as.numeric(details[model_nums[length(model_nums)],"df"])
    
    GOF = (1-pchisq(GLR, df=df_diff))
    return(GOF)
  }
  
  need_to_divide.AIC = function(model_nums,details,num_missing_models){
    # if models are missing (because y was the same for all the subset), logLik=0 and df=1 for
    # every missing model.
    sep_logLik = 0
    sep_df = num_missing_models
    model_nums = as.numeric(model_nums)
    for (r in model_nums[1:(length(model_nums)-1)]){
      sep_logLik = sep_logLik + as.numeric(details[r,"logLik"])
      sep_df = sep_df + as.numeric(details[r,"df"])
    }
    AIC.sep = 2*sep_df-2*sep_logLik
    AIC.all = as.numeric(details[model_nums[length(model_nums)],"AIC"])
    need_to_divide = AIC.sep<AIC.all
    # print(paste("AIC.sep = ",AIC.sep," AIC.all = ",AIC.all))
    return(need_to_divide)
  }
  
  get_theta = function(model.nb){
    tmp = strsplit(model.nb$family[[1]],"")
    tmp2 = tmp[[1]][19:(length(tmp[[1]])-1)]
    theta = as.numeric(paste(tmp2,collapse = ""))
    return(theta)
  }
  
  binary_matrix = function(data,input_iterate){
    tmp1 = NULL
    for (i in 1:length(input_iterate)){
      curr_col = input_iterate[i]
      vals = sort(unique(data[,curr_col]))
      print(curr_col)
      print(vals)
      tmp = matrix(0,dim(data)[1],(length(vals)+1))
      colnames(tmp) = paste(curr_col,1:(length(vals)+1),sep = ".")
      for (j in 1:length(vals)){
        tmp[,j] = data[,curr_col]==vals[j]
      }
      tmp[,(j+1)] = 1
      tmp1 = cbind(tmp1,tmp)
    }
    syn_exposure. = data[,"syn_exposure"]>0
    non_syn_exposure. = data[,"non_syn_exposure"]>0
    tmp1 = cbind(tmp1,syn_exposure.,non_syn_exposure.)
    return(as.matrix(tmp1))
  }
  
  process_data = function(curr_data,selected_input,selected_output,curr_details,count,output_num){
    curr_details[,colnames(selected_input)] = selected_input
    curr_details[,"output"] = output_num
    curr_details[,"rows_num"] = dim(curr_data)[1]
    curr_details[,"num_non_0_rows_1"] = (length(which(curr_data[,dim(curr_data)[2]]!=0))+1)
    curr_details[,"output_sum"] = sum(curr_data[,dim(curr_data)[2]])
    if (curr_details[,"output_sum"]==0 | dim(curr_data)[1]==1){
      curr_details[,c("logLik","P.logLik")] = 0
      curr_details[,c("df","P.df")] = c(2,1)
      curr_details[,c("AIC","P.AIC")] = c(4,2)
      curr_details[,c("theta","NB.converged","P.converged","P.test","P.GLR")] = NA
      curr_details[,"cols_num"] = dim(curr_data)[2]
    }else{   
      tmp = prepare_data_rm_cols(curr_data)
      curr_data = tmp[[1]]; curr_model_names = tmp[[2]]
      model.nb = find_model.nb.data(curr_data,curr_model_names)
      model.p = find_model.P.data(curr_data, curr_model_names)
      curr_details[,"cols_num"] = dim(curr_data)[2]
      if (is.double(model.nb)==FALSE){
        curr_details[,"logLik"] = logLik(model.nb)
        curr_details[,"df"] = find_df(model.nb) + 1 #+1 because theta is found outside of the glm
        curr_details[,"AIC"] = AIC(model.nb) + 2 #+2 because theta is found outside of the glm
        curr_details[,"theta"] = get_theta(model.nb)
        curr_details[,"NB.converged"] = model.nb$converged
      }else{
        print(selected_input)
        curr_details[,"logLik"] = NA
        curr_details[,"df"] = NA
        curr_details[,"AIC"] = NA
        curr_details[,"theta"] = NA
        curr_details[,"NB.converged"] = NA
      }
      if (is.double(model.p)==FALSE){
        if (is.double(model.nb)==FALSE){
          GLR.model.p = 2*(logLik(model.nb)[1]-logLik(model.p)[1])
          test.model.p = (1-pchisq(GLR.model.p, df=1))
        }else{
          GLR.model.p = NA
          test.model.p = NA
        }
        logLik.model.p = logLik(model.p)
        df.model.p = find_df(model.p)
        AIC.model.p = AIC(model.p)
        converged.p = model.p$converged 
      }else{
        GLR.model.p = NA
        test.model.p = NA
        logLik.model.p = NA 
        df.model.p = NA
        AIC.model.p = NA
        converged.p = NA
      }
      curr_details[,"P.GLR"] = GLR.model.p
      curr_details[,"P.test"] = test.model.p
      curr_details[,"P.logLik"] = logLik.model.p
      curr_details[,"P.df"] = df.model.p
      curr_details[,"P.AIC"] = AIC.model.p
      curr_details[,"P.converged"] = converged.p
    }
    return(curr_details)
  }
  
  is_neighbor_same_as_codon = function(curr_data,rm_cols,count){
    if (is.na(selected_input[1,"codon.pos"])==FALSE){
      if (selected_input[1,"codon.pos"]==1){
        rm_cols = c(rm_cols,which(colnames(curr_data)=="R.neighbor"))
      }
      if (selected_input[1,"codon.pos"]==2){
        rm_cols = c(rm_cols, 
                    which(colnames(curr_data)=="R.neighbor"),
                    which(colnames(curr_data)=="L.neighbor"))
      }
      if (selected_input[1,"codon.pos"]==3){
        rm_cols = c(rm_cols, which(colnames(curr_data)=="L.neighbor"))
      }
    }
    rm_cols = sort(unique(rm_cols))
    return(rm_cols)
  }
  
  loop_process = function(curr_data,curr_binary_data,selected_input,iterate_col,x){
    if (iterate_col!="syn_exposure" & iterate_col!="non_syn_exposure"){
      selected_input[,iterate_col] = x
    }
    curr_col = paste(iterate_col,x,sep = ".")
    if (is.null(dim(curr_binary_data)[1])==FALSE){
      plcs = which(curr_binary_data[,curr_col]==1)
      if (length(plcs)>1){
        curr_data = curr_data[plcs,]
        curr_binary_data = curr_binary_data[plcs,]
      }else{
        curr_data = NA
        curr_binary_data = NA
      }
    }else{
      curr_data = NA
      curr_binary_data = NA
    }
    tmp = list(selected_input,curr_binary_data,curr_data)
    return(tmp)
  }
  
  #############
  # DataSets: #
  #############
  output = c("syn","non_syn","transitions","transversions","y","A","C","G","T")
  
  # data:
  input = c("L.neighbor","R.neighbor","codon","site",
            "codon.pos","gene","mat_peptide","stem_loop","CG",             
            "amino_acid","exposure", "syn_exposure","non_syn_exposure")
  # input = c("site",
  #           # "codon",
  #           "regions","regions2",
  #           "codon.pos","backwords",
  #           # "CG",
  #           "CG_position",
  #           # "codon.a","codon.b","codon.c",
  #           # "R.neighbor","L.neighbor","amino_acid",
  #           "amino_acid_backwords",
  #           # "codon.a_backwords","codon.b_backwords","codon.c_backwords",
  #           "codon_backwords",
  #           "R.neighbor_backwords","L.neighbor_backwords",
  #           "exposure","syn_exposure","non_syn_exposure")
  # input_iterate = c("codon.pos","backwords","CG_position","codon_backwords","amino_acid_backwords",
  #                   "R.neighbor_backwords","L.neighbor_backwords","regions","regions2")
  
  input_iterate = c("L.neighbor","R.neighbor","codon","codon.pos","gene",
                    "mat_peptide","stem_loop","CG","amino_acid") 
  data = codons_table
  
  
  data_input = data[,input]
  data_output = data[,output]
  binary_data = binary_matrix(data,input_iterate)
  
  count = 0
  colnames_curr_details = c(input_iterate,"output","ID","logLik","df","AIC","theta","NB.converged","rows_num","cols_num",
                            "P.GLR","P.test","P.logLik","P.df","P.AIC","P.converged","num_non_0_rows_1","output_sum")
  details = matrix(0,500000,length(colnames_curr_details))
  selected_input = matrix(0,1,length(input_iterate))
  colnames(selected_input) = c(input_iterate)
  curr_details = matrix(0,1,dim(details)[2])
  colnames(details) = colnames(curr_details)
  colnames(curr_details) = colnames_curr_details
  iterate_vals = matrix(0,1,length(input_iterate))
  colnames(iterate_vals) = input_iterate
  for (i in 1:length(input_iterate)){
    curr_col = input_iterate[i]
    iterate_vals[,curr_col] = length(unique(data[,curr_col]))
  }
  iterate_vals = iterate_vals+2
  
  # save(iterate_vals,file = "iterate_vals2")
  
  parallel_vals = expand.grid(1:66,1:23,1:9)
  colnames(parallel_vals) = c("codon","amino_acid","output")
  q = parallel_vals[parallel_num,"codon"]
  w = parallel_vals[parallel_num,"amino_acid"]
  i = parallel_vals[parallel_num,"output"]
  
  
  # for (i in length(output)){
  selected_output = output[i]
  selected_input[1:length(selected_input)]=0
  curr_details[1:length(curr_details)]=0
  curr_data = cbind(data_input,data_output[,i])
  colnames(curr_data)[dim(curr_data)[2]] = selected_output
  rm_cols = NULL
  if (i<3){
    if (i==1){tmp = loop_process(curr_data,binary_data,selected_input,"syn_exposure","")}
    else{ #i=2
      tmp = loop_process(curr_data,binary_data,selected_input,"non_syn_exposure","")
    }
    selected_input = tmp[1][[1]]; curr_binary_data = tmp[2][[1]]; curr_data = tmp[3][[1]]
  }else{
    curr_binary_data = binary_data 
  }
  for (j in 1:iterate_vals[1,"codon.pos"]){
    curr_loop_data.j = curr_data; curr_loop_binary_data.j = curr_binary_data; rm_cols.j = rm_cols
    if (j==iterate_vals[1,"codon.pos"]){
      rm_cols = c(rm_cols,which(colnames(curr_data)=="codon.pos"))
      selected_input[,"codon.pos"] = NA
    }else{
      tmp = loop_process(curr_data,curr_binary_data,selected_input,"codon.pos",j)
      selected_input = tmp[1][[1]]; curr_binary_data = tmp[2][[1]]; curr_data = tmp[3][[1]]
      if(is.na(tmp[2])){
        curr_data = curr_loop_data.j; curr_binary_data = curr_loop_binary_data.j; rm_cols = rm_cols.j 
        next
      }
    }
    for (k in 1:iterate_vals[1,"mat_peptide"]){
      curr_loop_data.k = curr_data; curr_loop_binary_data.k = curr_binary_data; rm_cols.k = rm_cols
      if (k==iterate_vals[1,"mat_peptide"]){
        rm_cols = c(rm_cols,which(colnames(curr_data)=="mat_peptide"))
        selected_input[,"mat_peptide"] = NA
      }else{
        tmp = loop_process(curr_data,curr_binary_data,selected_input,"mat_peptide",k)
        selected_input = tmp[1][[1]]; curr_binary_data = tmp[2][[1]]; curr_data = tmp[3][[1]]
        if(is.na(tmp[2])){
          curr_data = curr_loop_data.k; curr_binary_data = curr_loop_binary_data.k; rm_cols = rm_cols.k 
          next
        }
      }
      for (p in 1:iterate_vals[1,"CG"]){
        curr_loop_data.p = curr_data; curr_loop_binary_data.p = curr_binary_data; rm_cols.p = rm_cols
        if (p==iterate_vals[1,"CG"]){
          rm_cols = c(rm_cols,which(colnames(curr_data)=="CG"))
          selected_input[,"CG"] = NA
        }else{
          tmp = loop_process(curr_data,curr_binary_data,selected_input,"CG",p)
          selected_input = tmp[1][[1]]; curr_binary_data = tmp[2][[1]]; curr_data = tmp[3][[1]]
          if(is.na(tmp[2])){
            curr_data = curr_loop_data.p; curr_binary_data = curr_loop_binary_data.p; rm_cols = rm_cols.p
            next
          }
        }
        for (a in 1:iterate_vals[1,"stem_loop"]){
          curr_loop_data.a = curr_data; curr_loop_binary_data.a = curr_binary_data; rm_cols.a = rm_cols
          if (a==iterate_vals[1,"stem_loop"]){
            rm_cols = c(rm_cols,which(colnames(curr_data)=="stem_loop"))
            selected_input[,"stem_loop"] = NA
          }else{
            tmp = loop_process(curr_data,curr_binary_data,selected_input,"stem_loop",a)
            selected_input = tmp[1][[1]]; curr_binary_data = tmp[2][[1]]; curr_data = tmp[3][[1]]
            if(is.na(tmp[2])){
              curr_data = curr_loop_data.a; curr_binary_data = curr_loop_binary_data.a; rm_cols = rm_cols.a
              next
            }
          }
          #for (q){#in 1:iterate_vals[1,"codon"]){
          curr_loop_data.q = curr_data; curr_loop_binary_data.q = curr_binary_data; rm_cols.q = rm_cols
          if (q==iterate_vals[1,"codon"]){
            rm_cols = c(rm_cols,which(colnames(curr_data)=="codon"))
            selected_input[,"codon"] = NA
            # use amino_acid instead of codons
            w_max = 1
          }else{
            tmp = loop_process(curr_data,curr_binary_data,selected_input,"codon",q)
            selected_input = tmp[1][[1]]; curr_binary_data = tmp[2][[1]]; curr_data = tmp[3][[1]]
            w_max = iterate_vals[1,"amino_acid"]
            if(is.na(tmp[2])){
              curr_data = curr_loop_data.q; curr_binary_data = curr_loop_binary_data.q; rm_cols = rm_cols.q
              next
            }
          }
          #for (w in w_max:iterate_vals[1,"amino_acid"]){
          w = max(w,w_max)
          curr_loop_data.w = curr_data; curr_loop_binary_data.w = curr_binary_data; rm_cols.w = rm_cols
          if (w==iterate_vals[1,"amino_acid"]){
            rm_cols = c(rm_cols,which(colnames(curr_data)=="amino_acid"))
            selected_input[,"amino_acid"] = NA
          }else{
            tmp = loop_process(curr_data,curr_binary_data,selected_input,"amino_acid",w)
            selected_input = tmp[1][[1]]; curr_binary_data = tmp[2][[1]]; curr_data = tmp[3][[1]]
            if(is.na(tmp[2])){
              curr_data = curr_loop_data.w; curr_binary_data = curr_loop_binary_data.w; rm_cols = rm_cols.w
              next
            }
          }
          if (j==1 |j==2){
            t_max=iterate_vals[1,"R.neighbor"]
          }else{t_max = 1}
          for (t in t_max:iterate_vals[1,"R.neighbor"]){
            curr_loop_data.t = curr_data; curr_loop_binary_data.t = curr_binary_data; rm_cols.t = rm_cols
            if (t == iterate_vals[1,"R.neighbor"]){
              rm_cols = c(rm_cols,which(colnames(curr_data)=="R.neighbor"))
              selected_input[,"R.neighbor"] = NA
            }else{
              tmp = loop_process(curr_data,curr_binary_data,selected_input,"R.neighbor",t)
              selected_input = tmp[1][[1]]; curr_binary_data = tmp[2][[1]]; curr_data = tmp[3][[1]]
              if(is.na(tmp[2])){
                curr_data = curr_loop_data.t; curr_binary_data = curr_loop_binary_data.t; rm_cols = rm_cols.t
                next
              }
            }
            if (j==2 |j==3){
              u_max=iterate_vals[1,"L.neighbor"]
            }else{u_max = 1}
            for (u in u_max:iterate_vals[1,"L.neighbor"]){
              curr_loop_data.u = curr_data; curr_loop_binary_data.u = curr_binary_data; rm_cols.u = rm_cols
              if (u == iterate_vals[1,"L.neighbor"]){
                rm_cols = c(rm_cols,which(colnames(curr_data)=="L.neighbor"))
                selected_input[,"L.neighbor"] = NA
              }else{
                tmp = loop_process(curr_data,curr_binary_data,selected_input,"L.neighbor",u)
                selected_input = tmp[1][[1]]; curr_binary_data = tmp[2][[1]]; curr_data = tmp[3][[1]]
                if(is.na(tmp[2])){
                  curr_data = curr_loop_data.u; curr_binary_data = curr_loop_binary_data.u; rm_cols = rm_cols.u
                  next
                }
              }
              for (m in 1:iterate_vals[1,"gene"]){
                curr_loop_data.m = curr_data; curr_loop_binary_data.m = curr_binary_data; rm_cols.m = rm_cols
                if (m==iterate_vals[1,"gene"]){
                  rm_cols = c(rm_cols,which(colnames(curr_data)=="gene"))
                  selected_input[,"gene"] = NA
                  # use regions2 instead of regions
                  n_max = 1
                }else{
                  tmp = loop_process(curr_data,curr_binary_data,selected_input,"gene",m)
                  selected_input = tmp[1][[1]]; curr_binary_data = tmp[2][[1]]; curr_data = tmp[3][[1]]; 
                  #n_max = iterate_vals[1,"regions2"]
                  if(is.na(tmp[2])){
                    curr_data = curr_loop_data.m; curr_binary_data = curr_loop_binary_data.m; rm_cols = rm_cols.m
                    next
                  }
                }
                # for (n in n_max:iterate_vals[1,"regions2"]){
                #   curr_loop_data.n = curr_data; curr_loop_binary_data.n = curr_binary_data; rm_cols.n = rm_cols
                #   if (n==iterate_vals[1,"regions2"]){
                #     rm_cols = c(rm_cols,which(colnames(curr_data)=="regions2"))
                #     selected_input[,"regions2"] = NA
                #   }else{
                #     tmp = loop_process(curr_data,curr_binary_data,selected_input,"regions2",n)
                #     selected_input = tmp[1][[1]]; curr_binary_data = tmp[2][[1]]; curr_data = tmp[3][[1]]
                #     if(is.na(tmp[2])){
                #       curr_data = curr_loop_data.n; curr_binary_data = curr_loop_binary_data.n; rm_cols = rm_cols.n
                #       next
                #     }
                #   }
                if (dim(curr_data)[1]!=0){
                  if(length(unique(curr_data[,dim(curr_data)[2]]))>1){
                    tmp = 1
                  }
                  rm_cols = c(rm_cols,which(colnames(curr_data)=="site"))
                  rel_cols = c(1:(dim(curr_data)[2]-4))
                  for (v in rel_cols){ # so that exposure/syn_exposure/non_syn_exposure will not be included!
                    if (length(unique(curr_data[,v]))==1){
                      rm_cols = c(rm_cols,v)
                    }
                  }
                  rm_cols = is_neighbor_same_as_codon(curr_data,rm_cols,count)
                  if (is.na(rm_cols)[1]==FALSE){
                    curr_data = curr_data[,-rm_cols]
                  }
                  if (is.na(dim(curr_data)[1])==FALSE){
                    count = count +1
                    print(count)
                    
                    curr_details = process_data(curr_data,selected_input,
                                                selected_output,curr_details,count,i)
                    print(curr_details)
                    details[count,] = as.numeric(curr_details)  
                  }
                  # }
                  # curr_data = curr_loop_data.n; curr_binary_data = curr_loop_binary_data.n; rm_cols = rm_cols.n
                }
                curr_data = curr_loop_data.m; curr_binary_data = curr_loop_binary_data.m; rm_cols = rm_cols.m
              }
              curr_data = curr_loop_data.u; curr_binary_data = curr_loop_binary_data.u; rm_cols = rm_cols.u
            }
            curr_data = curr_loop_data.t; curr_binary_data = curr_loop_binary_data.t; rm_cols = rm_cols.t
            # } 
            curr_data = curr_loop_data.w; curr_binary_data = curr_loop_binary_data.w; rm_cols = rm_cols.w
            # }
            curr_data = curr_loop_data.q; curr_binary_data = curr_loop_binary_data.q; rm_cols = rm_cols.q
          }
          curr_data = curr_loop_data.a; curr_binary_data = curr_loop_binary_data.a; rm_cols = rm_cols.a
        }
        curr_data = curr_loop_data.p; curr_binary_data = curr_loop_binary_data.p; rm_cols = rm_cols.p
      }
      curr_data = curr_loop_data.k; curr_binary_data = curr_loop_binary_data.k; rm_cols = rm_cols.k
    }   
    curr_data = curr_loop_data.j; curr_binary_data = curr_loop_binary_data.j; rm_cols = rm_cols.j
  }
  # }
  
  details = details[1:(count),]
  colnames(details) = colnames(curr_details)
  details = data.frame(details)
  file_name = paste0("details.q.",q,"_w.",w,"_i.",i,"_parallel_num.",parallel_num,collapse = "")
  setwd("/home/ubuntu/Dropbox/EC2_files/trees")
  save(details,file = file_name)
}

# print(Sys.time())
# details3 = create_details(2)
# print(Sys.time())

require(snow)
library(snow)
library(doSNOW)
no_cores =64
cl = makeCluster(no_cores,type="SOCK",outfile="")
registerDoSNOW(cl)
# codon_vals = 37:66
parallel_nums = c(1:13662)
print(Sys.time())
details = clusterApply(cl, parallel_nums,create_details)
print(Sys.time())
stopCluster(cl)


