############################
# Synonymous substitutions #
############################
# setwd("/Users/keren/Dropbox/covid/vars/")
# load("codons_table")
# transversion_to_transitions =
#   sum(codons_table$transversions)/sum(codons_table$transitions)
# save(transversion_to_transitions,file = "transversion_to_transitions")
# 
# edges = matrix("black",nrow = 64, ncol = 64)
# colnames(edges) <- c("AAA","AAC","AAG","AAT","ACA","ACC","ACG","ACT","AGA","AGC","AGG","AGT","ATA","ATC","ATG",
#                      "ATT","CAA","CAC","CAG","CAT","CCA","CCC","CCG","CCT","CGA","CGC","CGG","CGT","CTA","CTC",
#                      "CTG","CTT","GAA","GAC","GAG","GAT","GCA","GCC","GCG","GCT","GGA","GGC","GGG","GGT","GTA",
#                      "GTC","GTG","GTT","TAA","TAC","TAG","TAT","TCA","TCC","TCG","TCT","TGA","TGC","TGG","TGT",
#                      "TTA","TTC","TTG","TTT")
# rownames(edges) = colnames(edges)
# codon_names = colnames(edges)
# groups = edges
# g1 = c("GCT", "GCC", "GCA", "GCG")
# g2 = c("CGT", "CGC", "CGA", "CGG", "AGA", "AGG")
# g3 = c("AAT", "AAC")
# g4 = c("GAT", "GAC")
# g5 = c("TGT", "TGC")
# g6 = c("CAA", "CAG")
# g7 = c("GAA", "GAG")
# g8 = c("GGT", "GGC", "GGA", "GGG")
# g9 = c("CAT", "CAC")
# g10 = c("ATT", "ATC", "ATA")
# g11 = c("ATG")
# g12 = c("TTA", "TTG", "CTT", "CTC", "CTA", "CTG")
# g13 = c("AAA", "AAG")
# g14 = c("TTT", "TTC")
# g15 = c("CCT", "CCC", "CCA", "CCG")
# g16 = c("TCT", "TCC", "TCA", "TCG", "AGT", "AGC")
# g17 = c("ACT", "ACC", "ACA", "ACG")
# g18 = c("TGG")
# g19 = c("TAT", "TAC")
# g20 = c("GTT", "GTC", "GTA", "GTG")
# g21 = c("TAA", "TGA", "TAG")
# 
# # .1/2/3 number of possible mutations in a syn group in the i'th codon position
# group_sizes.3 = matrix(0,64,1)
# row.names(group_sizes.3) = codon_names
# group_sizes.2 = group_sizes.1 = group_sizes.3
# group_sizes.3[c(g1,g8,g9,g15,g16,g17,g20),1] = 1+2*ratio
# group_sizes.3[c(g2,g12),1] = 1+2*ratio
# group_sizes.3[c("AGA", "AGG","TAA","TAG"),1] = 1
# group_sizes.3[c(g3,g4,g5,g6,g7,g13,g14,g19),1] = 1
# group_sizes.3[c(g10),1] = 1+ratio
# group_sizes.3["ATA",1] = ratio
# group_sizes.3[c("TTA", "TTG","AGT", "AGC"),1] = 1
# 
# group_sizes.2[c("TAA", "TGA"),1] = 1
# group_sizes.1[c("CGA", "CGG", "AGA", "AGG"),1] = ratio
# group_sizes.1[c("TTA", "TTG", "CTA", "CTG"),1] = 1
# 
# 
# group_names = c("Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly","His","Ile","Start/Met",
#                 "Leu","Lys","Phe","Pro","Ser","Thr","Trp","Tyr","Val","Stop")
# group_sizes = c(4,6,2,2,2,2,2,4,2,3,1,6,2,2,4,6,4,1,2,4,3)
# 
# group_names2 = c("A","R","N","D","C","Q","E","G","H","I","M",
#                  "L","K","F","P","S","T","W","Y","V","Stop")
# codons_to_amino_acids = matrix(0,64,2)
# row.names(codons_to_amino_acids) = codon_names
# for (j in 1:length(group_names)){
#   gr = paste("g",j,sep="")
#   curr.g = eval(parse(text = gr))
#   codons_to_amino_acids[curr.g,1] = group_names[j]
#   codons_to_amino_acids[curr.g,2] = group_names2[j]
# }

# 
# 
# save(groups,file = "groups")
# save(group_sizes,file = "group_sizes")
# save(codons_to_amino_acids,file = "codons_to_amino_acids")
# save(group_sizes.1,file = "group_sizes.1")
# save(group_sizes.2,file = "group_sizes.2")
# save(group_sizes.3,file = "group_sizes.3")


##############
# Functions: #
##############
add_CG = function(data){
  CG = matrix(0,dim(data)[1],1)
  CG[which(data$base=="C" & data$R.neighbor=="G")] = 1
  CG[which(data$base=="G" & data$L.neighbor=="C")] = 2
  data = cbind(data,CG)
  return(data)
}

add_transitions_transversions = function(data){
  transitions = transversions = matrix(0,dim(data),1)
  transitions[which(data$base=="A")] = data$G[which(data$base=="A")]
  transitions[which(data$base=="G")] = data$A[which(data$base=="G")]
  transitions[which(data$base=="C")] = data$T[which(data$base=="C")]
  transitions[which(data$base=="T")] = data$C[which(data$base=="T")]
  
  transversions[which(data$base=="A")] = data$C[which(data$base=="A")]+
    data$T[which(data$base=="A")]
  transversions[which(data$base=="C")] = data$A[which(data$base=="C")]+
    data$G[which(data$base=="C")]
  transversions[which(data$base=="G")] = data$C[which(data$base=="G")]+
    data$T[which(data$base=="G")]
  transversions[which(data$base=="T")] = data$A[which(data$base=="T")]+
    data$G[which(data$base=="T")]
  
  return(cbind(data,transitions,transversions))
}


add_amino_acid_exposure = function(data,codons_to_amino_acids,ratio,group_sizes.1,
                                   group_sizes.2,group_sizes.3){
  amino_acid = syn_exposure =  non_syn_exposure = matrix(0,dim(data)[1],1)
  # plcs = which(is.element(data$codon,rownames(codons_to_amino_acids)))
  # no_valid_codon_plcs = which(is.element(data$codon,rownames(codons_to_amino_acids)==FALSE))
  for (i in 1:dim(codons_to_amino_acids)[1]){
    print(i)
    curr_codon = rownames(codons_to_amino_acids)[i]
    amino_acid[which(data$codon==curr_codon)] = codons_to_amino_acids[i,2]
    syn_exposure[which(data$codon==curr_codon & data$codon.pos==1)] = 
      group_sizes.1[curr_codon,1]
    syn_exposure[which(data$codon==curr_codon & data$codon.pos==2)] = 
      group_sizes.2[curr_codon,1]
    syn_exposure[which(data$codon==curr_codon & data$codon.pos==3)] = 
      group_sizes.3[curr_codon,1]
  }
  non_syn_exposure = (1+2*ratio)-syn_exposure
  return(cbind(data,amino_acid,syn_exposure,non_syn_exposure))
}

add_syn_non_syn = function(data){
  output_codons = matrix(0,dim(data)[1],4)
  input_codons = unlist(data$codon)
  colnames(output_codons) = c("A","C","G","T")
  output_codons = data.frame(output_codons)
  data$codon.pos = as.numeric(data$codon.pos)
  for (i in 1:dim(output_codons)[1]){
    #print(i)
    for(j in 1:4){
      tmp = strsplit(input_codons[i],split = "")[[1]]
      tmp[data$codon.pos[i][[1]]] = colnames(output_codons)[j]
      output_codons[i,j] = paste0(tmp,collapse = "")
    }
  }
  output_amino_acid = matrix(0,dim(data)[1],4)
  colnames(output_amino_acid) = c("A.amino","C.amino","G.amino","T.amino")
  syn_flag = output_amino_acid
  plcs = which(is.element(input_codons,row.names(codons_to_amino_acids)))
  for (i in 1:4){
    output_amino_acid[plcs,i] = codons_to_amino_acids[output_codons[plcs,i],2]
    syn_flag[which(data$amino_acid==output_amino_acid[,i] & data$amino_acid!=0),i] = 1
  }
  syn = apply(syn_flag*data[,c("A","C","G","T")],1,sum)
  non_syn = apply( (syn_flag-1)*(-1) * data[,c("A","C","G","T")],1,sum)
  syn[-plcs] = NA
  non_syn[-plcs] = NA
  return(cbind(data,syn,non_syn,output_amino_acid))
}

unlist_table = function(old_table){
  new_table = NULL
  for (i in c(1:dim(old_table)[2])){
    print(i)
    tmp = unlist(old_table[,i])
    new_table = cbind(new_table,tmp)
  }
  colnames(new_table) = colnames(old_table)
  new_table = data.frame(new_table)
  return(new_table)
}

change_cols_to_numeric = function(numeric_cols,data){
  for (i in 1:length(numeric_cols)){
    data[,numeric_cols[i]] = as.numeric(data[,numeric_cols[i]])
  }
  return(data)
}

##############
setwd("/Users/keren/Dropbox/covid/vars/tables_data/")
load("base_outputs")
load("base_states")
load("codon_states")
load("codon_outputs")

setwd("/Users/keren/Dropbox/covid/vars/")
load("site_details")
load("groups")
load("group_sizes")
load("codons_to_amino_acids")
load("group_sizes.1")
load("group_sizes.2")
load("group_sizes.3")
load("transversion_to_transitions")
ratio = transversion_to_transitions

setwd("/Users/keren/Dropbox/EC2_files/")
source("COVID_functions.R")
tree = read.tree("tree.nwk")

tmp = strsplit(codon_states,split = "")
codons_table = NULL
for (i in 1:length(tmp)){
  print(i)
  if (tmp[i][[1]][7]=="s"){ # codon is "missing"
    site = as.numeric(paste0(tmp[i][[1]][11:length(tmp[i][[1]])],collapse = ""))
  }else{
    site = as.numeric(paste0(tmp[i][[1]][7:length(tmp[i][[1]])],collapse = ""))
  }
  # new_line = c(tmp[i][[1]][1],tmp[i][[1]][2],tmp[i][[1]][3],
  #              paste0(tmp[i][[1]][4:6],collapse = ""),
  #              site)
  # new_line = c(new_line,regions[regions$site==site,])
  # regions_table = rbind(regions_table,regions[regions$site==site,])
  codons_table = rbind(codons_table,
                       c(tmp[i][[1]][1],tmp[i][[1]][2],tmp[i][[1]][3],
                         paste0(tmp[i][[1]][4:6],collapse = ""),
                         site,regions[regions$site==site,]))
}
missing_codons_plcs = which(codons_table$codon=="mis")
colnames(codons_table)[1:5] = c("L.neighbor","R.neighbor","base","codon","site")
codons_table = cbind(codons_table,codon_outputs[1:dim(codons_table)[1],])
codons_table = data.frame(codons_table)
# codons_table = codons_table[which(codons_table$L.neighbor!="?" & 
#                                     codons_table$R.neighbor!="?"),]
codons_table = add_CG(codons_table)
codons_table = add_transitions_transversions(codons_table)
codons_table = unlist_table(codons_table)
codons_table = add_amino_acid_exposure(codons_table,
                                       codons_to_amino_acids,ratio,group_sizes.1,group_sizes.2,group_sizes.3)

numeric_cols = c("site","site_pattern","ref_site","codon.pos",
                 "A","C","G","T","line","exposure","zero_branch",
                 "CG","transitions","transversions","syn_exposure",
                 "non_syn_exposure")
codons_table = change_cols_to_numeric(numeric_cols,data = codons_table)

rm_cols = which(is.element(colnames(codons_table),c("M","q","s","site.1","UTR")))
codons_table = codons_table[,-rm_cols]
y = codons_table$transitions+codons_table$transversions
codons_table = cbind(codons_table,y)
codons_table = add_syn_non_syn(codons_table)

setwd("/Users/keren/Dropbox/covid/vars/")
save(codons_table,file = "codons_table2")
load("codons_table2")

# #check:
# plcs = NULL
# tmp = strsplit(codons_table$codon,split = "")
# for (i in 1: dim(codons_table)[1]){
#   if(tmp[[i]][codons_table$codon.pos[i]]!=codons_table$base[i]){
#     plcs = c(plcs,i)
#   }
# }

tmp = strsplit(base_states,split = "")
base_table = NULL
for (i in 1:length(tmp)){
  print(i)
  site = as.numeric(paste0(tmp[i][[1]][4:length(tmp[i][[1]])],collapse = ""))
  new_line = c(tmp[i][[1]][1],tmp[i][[1]][2],tmp[i][[1]][3],site)
  new_line = c(new_line,regions[regions$site==site,])
  base_table = rbind(base_table,new_line)
}
colnames(base_table)[1:4] = c("L.neighbor","R.neighbor","base","site")
base_table = cbind(base_table,base_outputs[1:dim(base_table)[1],])
base_table = data.frame(base_table)

base_table = base_table[which(base_table$base!="?"),]
base_table = add_CG(base_table)
base_table = add_transitions_transversions(base_table)
rm_cols = which(is.element(colnames(base_table),c("M","q","s","site.1")))
base_table = base_table[,-rm_cols]
base_table = unlist_table(base_table)
numeric_cols = c("site","site_pattern","ref_site",
                 "A","C","G","T","line","exposure","zero_branch",
                 "CG","transitions","transversions")
base_table = change_cols_to_numeric(numeric_cols,data = base_table)
y = base_table$transitions+base_table$transversions
base_table = cbind(base_table,y)

setwd("/Users/keren/Dropbox/covid/vars/")
save(base_table,file = "base_table2")
load("base_table2")

plcs = NULL
for (i in 1:dim(codons_table)[1]){
  tmp = strsplit(codons_table$codon[i],split = "")[[1]]
  if (tmp[codons_table$codon.pos[i]]!=codons_table$base[i]){
    plcs = c(plcs,i)
  }
}
