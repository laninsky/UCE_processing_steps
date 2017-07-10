# Concatenation of multiple fasta files. Spits out a concatenated fasta and phylip file, and a tab-delimited partition file listing the 
# position of these partitions. Requires the path and name of a text file which has the paths and file names of the fasta files you are 
# interested in concatenating

# change the file name read in for the read.table function
file_names <- as.matrix(read.table("C:\\Users\\Alana\\Dropbox\\beetles\\mitogenomes_18S_28S\\18S_28S\\aligned\\file.txt"))

no_of_files <- dim(file_names)[1]

intable <- read.table(file_names[1,1],header=FALSE,stringsAsFactors=FALSE,sep="\t")
rows <- dim(intable)[1]
to_write_final <- intable[1,1]
sequencepaste <- NULL
for (j in 2:rows) {
  if ((length(grep(">",intable[j,1])))>0) {
    to_write <- toupper(sequencepaste)
    to_write_final <- rbind(to_write_final,to_write)
    to_write_final <- rbind(to_write_final,intable[j,1])
    sequencepaste <- NULL
  } else {
    sequencepaste <- paste(sequencepaste,intable[j,1],sep="")
  }
}
final_fasta <- rbind(to_write_final,toupper(sequencepaste))
final_partition <- c(file_names[1,1],1,nchar(final_fasta[2,1]))
final_partition <- t(as.matrix(final_partition))


for (i in 2:no_of_files) {
  intable <- read.table(file_names[i,1],header=FALSE,stringsAsFactors=FALSE,sep="\t")
  rows <- dim(intable)[1]
  to_write_final <- intable[1,1]
  sequencepaste <- NULL
    
  for (j in 2:rows) {
    if ((length(grep(">",intable[j,1])))>0) {
      to_write <- toupper(sequencepaste)
      to_write_final <- rbind(to_write_final,to_write)
      to_write_final <- rbind(to_write_final,intable[j,1])
      sequencepaste <- NULL
    } else {
      sequencepaste <- paste(sequencepaste,intable[j,1],sep="")
    }
  }

  to_write_final <- rbind(to_write_final,toupper(sequencepaste))
  final_partition <- rbind(final_partition,c(file_names[i,1],(as.numeric(final_partition[(i-1),3])+1),(as.numeric(final_partition[(i-1),3])+nchar(to_write_final[2,1]))))
  nchar_so_far <- nchar(final_fasta[2,1])
  
  for (j in seq(1,dim(to_write_final)[1],2)) {
    if(to_write_final[j,1] %in% final_fasta[(seq(1,dim(final_fasta)[1],2)),1]) {
      k <- which(final_fasta[,1]==to_write_final[j,1])
      final_fasta[(k+1),1] <- paste(final_fasta[(k+1),1],to_write_final[(j+1),1],sep="")
    } else {
      final_fasta <- rbind(final_fasta,to_write_final[j,1])
      tempseq <- paste(rep("-",nchar_so_far),collapse="")
      final_fasta <- rbind(final_fasta,(paste(tempseq,to_write_final[(j+1),1],sep="")))
    }
  }  
  
  for (j in seq(2,dim(final_fasta)[1],2)) {
    if(nchar(final_fasta[j,1])<as.numeric(final_partition[(dim(final_partition)[1]),3])) {
      no_new_chars <- nchar(to_write_final[2,1])
      tempseq <- paste(rep("-",no_new_chars),collapse="")
      final_fasta[j,1] <- paste(final_fasta[j,1],tempseq,sep="")
    }
  }
}

phylip <- cbind(final_fasta[(seq(1,(dim(final_fasta)[1]),2)),1],final_fasta[(seq(2,(dim(final_fasta)[1]),2)),1])
phylip[,1] <- gsub(">","",phylip[,1],fixed=TRUE)
first_phylip_line <- c(dim(phylip)[1],final_partition[(dim(final_partition)[1]),3])
phylip <- rbind(first_phylip_line,phylip)

write.table(final_fasta,"concatenated_fasta.fasta",quote=FALSE, row.names=FALSE,col.names=FALSE)
write.table(phylip,"concatenated_phylip.phylip",quote=FALSE, row.names=FALSE,col.names=FALSE)
write.table(final_partition,"partition_file.txt",quote=FALSE, row.names=FALSE,col.names=FALSE)



