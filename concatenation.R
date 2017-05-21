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
    } else {
      sequencepaste <- paste(sequencepaste,intable[j,1],sep="")
    }
  }

  to_write_final <- rbind(to_write_final,toupper(sequencepaste))
  
  for (j in seq(1,dim(to_write_final)[1],2) {

}
  
