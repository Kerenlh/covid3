# setwd("C:/covid")
setwd("/Users/keren/Dropbox/covid")
source("COVID_functions.R")

tree = read.tree("tree.nwk")
tree = short_names_tree(tree)

# accession_numbers = NULL
setwd("/Users/keren/Dropbox/covid/tree_sequences")
seqs = NULL
for (i in 259:length(tree$tip.label)){
  print(i)
  # tmp = strsplit( tree$tip.label[[i]],split = "\\|")
  # accession_numbers = c(accession_numbers,tmp[[1]][1])
  result = tryCatch({
    # seqs = c(seqs,read.GenBank(tree$tip.label[i], seq.names = tree$tip.label[i], species.names = TRUE,
    #                            as.character = TRUE, chunk.size = 400, quiet = TRUE))
    seqs = read.GenBank(tree$tip.label[i], seq.names = tree$tip.label[i], species.names = TRUE,
                               as.character = TRUE, chunk.size = 400, quiet = TRUE)
    cat(file=paste0(tree$tip.label[i],".fasta",collapse = ""), paste(paste0(">",names(seqs)),
                                     sapply(seqs, paste, collapse=""), sep="\n"), sep="\n");
    rm(seqs)
  }, warning = function(w) {
  }, error = function(e) {
    error_seqs = c(error_seqs,i)
    print(paste("error in: ",i))
  }, finally = {
  })
}
  
  
  # accession_number = tmp[[1]][1]
  # file_name = paste0(accession_number,".fasta")
  # seq = read.GenBank(accession_number, seq.names = accession_number, species.names = TRUE,
  #                            as.character = TRUE, chunk.size = 400, quiet = TRUE)
  # cat(file=file_name, paste(paste0(">",names(seq)),
  #                              sapply(seqs, paste, collapse=""), sep="\n"), sep="\n");


# cat(file="all_seqs.fasta", paste(paste0(">",names(seqs)),
#                              sapply(seqs, paste, collapse=""), sep="\n"), sep="\n");
# seqs_fatsa = readAAStringSet("all_seqs.fasta")

# tmp = readAAStringSet(file_name)
# tmp = readDNAStringSet(file_name)
aligned_seqs = msa(seqs_fatsa)
# tmp = as.DNAbin.list(seqs)
# aligned_seqs = msa(seqs)