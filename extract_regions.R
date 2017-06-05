# 

working_dir <- "C:\\Users\\Alana\\Dropbox\\beetles\\mitogenomes_18S_28S\\mitogenome\\final"
fasta_suffix <- ".fasta"
feature_suffix <- "_feature_table.txt"

setwd(working_dir)

for (i in list.files(pattern=fasta_suffix)) {
    temp_fasta <-  read.table(i,header=FALSE,stringsAsFactors=FALSE,sep="\t")
    fasta <- temp_fasta[1,1]
    sequencepaste <- NULL
    for (j in 2:(dim(temp_fasta)[1])) {
      if ((length(grep(">",temp_fasta[j,1])))>0) {
      to_write <- toupper(sequencepaste)
      fasta <- rbind(fasta,to_write)
      fasta <- rbind(fasta,temp_fasta[j,1])
      sequencepaste <- NULL
      } else {
      sequencepaste <- paste(sequencepaste,temp_fasta[j,1],sep="")
      }
    }
    fasta <- rbind(fasta,toupper(sequencepaste))
