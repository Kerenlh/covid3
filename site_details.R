setwd("/Users/keren/Dropbox/covid")
source("COVID_functions.R")

setwd("/Users/keren/Dropbox/EC2_files/")
aligned_seqs = readDNAStringSet(file="aligned_tree_seqs_decipher")
tree = read.tree("tree.nwk")

ref_seq = "NC_045512.2"
plcs = which(strsplit(as.character(aligned_seqs[[ref_seq]]),split = "")[[1]]!="-")
ref_regions = matrix(0,length(plcs),8)
colnames(ref_regions) = c("ref_site","codon.pos","UTR","gene",
                          "mat_peptide","stem_loop","notes","notes2")
ref_regions = data.frame(ref_regions)
ref_regions$ref_site = 1:length(plcs)


ref_regions$UTR[1:265] = "5UTR"
ref_regions$gene[266:21555] = "ORF1ab"
ref_regions$codon.pos[13468:21555] = matrix(1:3, length(13468:21555),1)
ref_regions$codon.pos[266:13468] = matrix(1:3, length(266:13468),1)
ref_regions$mat_peptide[266:805] = "ORF1ab"
ref_regions$notes[266:805] = "nsp1"
ref_regions$mat_peptide[806:2719] = "ORF1ab"
ref_regions$notes[806:2719] = "nsp2"
ref_regions$notes[2720:8554] = "nsp3"
ref_regions$mat_peptide[2720:8554] = "ORF1ab"
ref_regions$notes[8555:10054] = "nsp4"
ref_regions$mat_peptide[8555:10054] = "ORF1ab"
ref_regions$notes[10055:10972] = "nsp5"
ref_regions$mat_peptide[10055:10972] = "ORF1ab"
ref_regions$notes[10973:11842] = "nsp6"
ref_regions$mat_peptide[10973:11842] = "ORF1ab"
ref_regions$notes[11843:12091] = "nsp7"
ref_regions$mat_peptide[11843:12091] = "ORF1ab"
ref_regions$notes[12092:12685] = "nsp8"
ref_regions$mat_peptide[12092:12685] = "ORF1ab"
ref_regions$notes[12686:13024] = "nsp9"
ref_regions$mat_peptide[12686:13024] = "ORF1ab"
ref_regions$notes[13025:13441] = "nsp10"
ref_regions$mat_peptide[13025:13441] = "ORF1ab"
ref_regions$notes[13442:16236] = "nsp12"
ref_regions$mat_peptide[13442:16236] = "ORF1ab"
ref_regions$notes[16237:18039] = "nsp13"
ref_regions$mat_peptide[16237:18039] = "ORF1ab"
ref_regions$notes[18040:19620] = "nsp14"
ref_regions$mat_peptide[18040:19620] = "ORF1ab"
ref_regions$notes[19621:20658] = "nsp15"
ref_regions$mat_peptide[19621:20658] = "ORF1ab"
ref_regions$notes[20659:21552] = "nsp16"
ref_regions$mat_peptide[20659:21552] = "ORF1ab"
ref_regions$codon.pos[266:13483] = matrix(1:3, length(266:13483),1)
ref_regions$notes2[13442:13480] = "nsp11"
ref_regions$mat_peptide[13442:13480] = "ORF1ab"
ref_regions$stem_loop[13476:13503] = "ORF1ab"
ref_regions$notes[13476:13503] = "stem_loop_1"
ref_regions$notes[13488:13542] = "stem_loop_2"
ref_regions$stem_loop[13488:13542] = "ORF1ab"
ref_regions$gene[21563:25384] = "S"
ref_regions$gene[25393:26220] = "ORF3a"
ref_regions$codon.pos[21563:25384] = matrix(1:3, length(21563:25384),1)
ref_regions$codon.pos[25393:26220] = matrix(1:3, length(25393:26220),1)
ref_regions$gene[26245:26472] = "E"
ref_regions$codon.pos[26245:26472] = matrix(1:3, length(26245:26472),1)
ref_regions$gene[26523:27191] = "M"
ref_regions$codon.pos[26523:27191] = matrix(1:3, length(26523:27191),1)
ref_regions$gene[27202:27387] = "ORF6"
ref_regions$codon.pos[27202:27387] = matrix(1:3, length(27202:27387),1)
ref_regions$gene[27394:27759] = "ORF7a"
ref_regions$codon.pos[27756:27887] = matrix(1:3, length(27756:27887),1)
ref_regions$codon.pos[27394:27759] = matrix(1:3, length(27394:27759),1)
ref_regions$gene[27756:27758] = "ORF7a_ORF7b"
ref_regions$gene[27756:27887] = "ORF7b"
ref_regions$gene[27894:28259] = "ORF8"
ref_regions$codon.pos[27894:28259] = matrix(1:3, length(27894:28259),1)
ref_regions$gene[28274:29533] = "N"
ref_regions$codon.pos[28274:29533] = matrix(1:3, length(28274:29533),1)
ref_regions$gene[29558:29674] = "ORF10"
ref_regions$codon.pos[29558:29674] = matrix(1:3, length(29558:29674),1)
ref_regions$stem_loop[29609:29644] = "ORF10"
ref_regions$notes[29609:29644] = "stem_loop_1"
ref_regions$stem_loop[29629:29657] = "ORF10"
ref_regions$notes2[29629:29657] = "stem_loop_2"
ref_regions$UTR[29675:29903] = "3UTR"
ref_regions$stem_loop[29728:29768] = "stem_loop_2_like"

regions = matrix(0,length(aligned_seqs[1][[1]]),9)
colnames(regions) = c("site",colnames(ref_regions)) 
regions = data.frame(regions)
regions[plcs,2:9]=ref_regions
regions$site = 1:dim(regions)[1]

setwd("/Users/keren/Dropbox/covid/vars")
load("all_site_patterns")
load("all_no_subs_plcs")
load("all_site_patterns_plcs")
# load("all_subs")

unique_site_patterns = unique(all_site_patterns)
site_pattern_plcs_indices = matrix(0,length(all_site_patterns_plcs),1)
for (i in 1:length(unique_site_patterns)){
#  print(i)
  plcs = which(all_site_patterns==unique_site_patterns[i])
  site_pattern_plcs_indices[plcs] = i
}

site_pattern = matrix(0,(length(all_site_patterns_plcs)+length(all_no_subs_plcs)),1)
site_pattern[all_site_patterns_plcs] = site_pattern_plcs_indices
which(site_pattern[all_no_subs_plcs]!=0) #check

load("missing_all_site_patterns")
load("missing_all_site_patterns_plcs")
load("missing_all_no_subs_plcs")
# all_site_patterns_plcs = unique(all_site_patterns_plcs)

unique_site_patterns = unique(all_site_patterns)
site_pattern_plcs_indices = matrix(0,length(all_site_patterns_plcs),1)
for (i in 1:length(unique_site_patterns)){
  #  print(i)
  plcs = which(all_site_patterns==unique_site_patterns[i])
  site_pattern_plcs_indices[plcs] = -i
}

# site_pattern = matrix(0,(length(all_site_patterns_plcs)+length(all_no_subs_plcs)),1)
site_pattern[all_site_patterns_plcs] = site_pattern_plcs_indices
which(site_pattern[all_no_subs_plcs]!=0) #check


regions = cbind(site_pattern,regions)
setwd("/Users/keren/Dropbox/covid/vars/")
save (regions,file = "site_details")

load("site_details")

ref_seq = read.GenBank("NC_045512.2")
tmp = NULL
for (i in 1:length(ref_seq$NC_045512.2)){
  if (ref_seq$NC_045512.2[i]=="88"){
    letter = "A"
  }else if (ref_seq$NC_045512.2[i]=="18"){
    letter = "T"
  }else if(ref_seq$NC_045512.2[i]=="48"){
    letter = "G"
  }else if(ref_seq$NC_045512.2[i]=="28"){
    letter = "C"
  }
  tmp = c(tmp,letter)
}

ref_seq = matrix(0,dim(regions)[1],1)
ref_seq[which(regions$ref_site>0)] = tmp
regions = cbind(regions,ref_seq)

codons = NULL
k=1
curr_codon_sites = which(regions$codon.pos>0)
curr_codon.pos = regions$codon.pos[curr_codon_sites]
seq = regions$ref_seq[unique(curr_codon_sites)]
while(k <=length(curr_codon_sites)){
  if (sum(curr_codon.pos[k:(k+2)]==c(1,2,3))==3){
    codons = c(codons,rep(paste0(seq[k:(k+2)],collapse = ""),3))
    k = k+3
  }else{
    codons = c(codons,"missing")
    print(paste("k=",k,"i=",i))
    k= k+1
  }
}

ref_codon = matrix(0,dim(regions)[1],1)
ref_codon[curr_codon_sites] = codons
regions = cbind(regions,ref_codon)

load("codons_to_amino_acids")
amino_acid = matrix(0,dim(regions)[1],1)
for (i in 1:dim(codons_to_amino_acids)[1]){
  print(i)
  curr_codon = rownames(codons_to_amino_acids)[i]
  amino_acid[which(regions$ref_codon==curr_codon)] = codons_to_amino_acids[i,2]
}

ref_amino_acid = amino_acid
regions = cbind(regions,ref_amino_acid)

setwd("/Users/keren/Dropbox/covid/vars/")
save (regions,file = "site_details")
load("site_details")

########
regions$codon.pos[11118:11120] = c(1,2,3)
regions$codon.pos[11123:11125] = c(3,1,2)
regions$codon.pos[11132:11134] = c(3,1,2)
regions$codon.pos[11137:11139] = c(2,3,1)

regions$gene[11118:11120] = "ORF1ab"
regions$gene[11123:11125] = "ORF1ab"
regions$gene[11132:11134] = "ORF1ab"
regions$gene[11137:11139] = "ORF1ab"

regions$mat_peptide[11118:11120] = "ORF1ab"
regions$mat_peptide[11123:11125] = "ORF1ab"
regions$mat_peptide[11132:11134] = "ORF1ab"
regions$mat_peptide[11137:11139] = "ORF1ab"
regions$notes[11118:11120] = "nsp6"
regions$notes[11123:11125] = "nsp6"
regions$notes[11132:11134] = "nsp6"
regions$notes[11137:11139] = "nsp6"


regions$codon.pos[22411:22419] = c(2,3,1,2,3,1,2,3,1)
regions$gene[22411:22419] = "S"

regions$codon.pos[25759:25761] = c(3,1,2)
regions$gene[25759:25761] = "ORF3a"

regions$codon.pos[21441:21443] = c(2,3,1)
regions$gene[21441:21443] = "ORF1ab"
regions$mat_peptide[21441:21443] = "ORF1ab"
regions$notes[21441:21443] = "nsp16"

regions$codon.pos[26156:26158] = c(1,2,3)
regions$gene[26156:26158] = "ORF3a"

regions$codon.pos[28324:28326] = c(1,2,3)
regions$gene[28324:28326] = "ORF8"

setwd("/Users/keren/Dropbox/covid/vars/")
save (regions,file = "site_details")
setwd("/Users/keren/Dropbox/EC2_files/")
save (regions,file = "site_details")



#gaps_in_coding = c(11118:11120,11123:11125,11132:11134,11137:11139,
                   # 21441:21443,22411:22419,25759:25761,26156:26158,
                   # 28324:28326)
# 29646:29647: This is also a gap of 2, an insertion appears only in one seq, 
# relevant node is i=5350,  parent_name = "no_name_14339", child_name "MT520337.1"


codon_pos_1_plcs = which(regions$codon.pos==1)
codon_pos_2_plcs = which(regions$codon.pos==2)
codon_pos_3_plcs = which(regions$codon.pos==3)
codon_pos_plcs = cbind(codon_pos_1_plcs[4411:length(codon_pos_1_plcs)],
                       codon_pos_2_plcs[4412:length(codon_pos_2_plcs)],
                       codon_pos_3_plcs[4412:length(codon_pos_3_plcs)])
colnames(codon_pos_plcs) = c("pos1","pos2","pos3")
codon_pos_plcs = data.frame(codon_pos_plcs)
codon_pos_plcs = cbind(codon_pos_plcs$pos1[4729:length(codon_pos_plcs$pos1)],
                       codon_pos_plcs$pos2[4730:length(codon_pos_plcs$pos2)],
                       codon_pos_plcs$pos3[4730:length(codon_pos_plcs$pos3)])
colnames(codon_pos_plcs) = c("pos1","pos2","pos3")
codon_pos_plcs = data.frame(codon_pos_plcs)

tmp = which(codon_pos_plcs$pos2-codon_pos_plcs$pos1!=1)
tmp = which(codon_pos_3_plcs-codon_pos_1_plcs!=2)

codon_plcs = which(regions$codon.pos!=0)

# Sites 13468:13484 have double codon meaning - it should be 13468:21555 
# but it overlaps with 266:13483 so 13468 is kept as codon.pos=3 but 
# according to 13468:21555 it should be 1.
# Sites 27756:27759 have the same problem (27394:27759,27756:27887)

###########
# Check : #
load("ORF1ab1")
plcs = which(regions$gene=="ORF1ab")
ncbi_amino_acid = matrix(0,dim(regions),1)
j=1
count = 0
for (i in 1:length(plcs)){
  count = count+1
  ncbi_amino_acid[plcs[i]] = ORF1ab1[j]
  if (count==3){
    count = 0
    j = j+1
  }
}

check = cbind(ncbi_amino_acid,regions$ref_amino_acid)
colnames(check) = c("ncbi","regions")
plcs = which(check[,1]!=check[,2] & check[,1]!="0")

load("ORF1ab2")
plcs = which(regions$gene=="ORF1ab")
j=1
count = 0
for (i in 12244:length(plcs)){
  count = count+1
  ncbi_amino_acid[plcs[i]] = ORF1ab2[j]
  if (count==3){
    count = 0
    j = j+1
  }
}

ncbi_amino_acid = ncbi_amino_acid[-13524]
check = cbind(ncbi_amino_acid,regions$ref_amino_acid,regions$ref_amino_acid2)
colnames(check) = c("ncbi","regions","codon.pos2")
plcs = which(check[,1]!=check[,2] & check[,1]!="0")
plcs = which(check[,1]!=check[,3] & check[,1]!="0")

ref_codon.pos2 = ref_regions$codon.pos
ref_codon.pos2[13471:21555] = matrix(1:3,length(13471:21555))
codon.pos2 = matrix(0,dim(regions)[1],1)
codon.pos2[which(regions$ref_site!="0")] = ref_codon.pos2
           
