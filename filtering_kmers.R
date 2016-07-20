#To run this script modify the path to "overrepresented" to your fasta file containing the overrepresented kmers across your samples
#If the kmer starts in the second position of the read etc prefix it with an 'N' e.g. 'NTGCCTC'
#The script should be placed in a folder with your F and R reads and executed by `Rscript filtering_kmers.R`
#In order for the script to work, you should have a R1 and R2 file for each individual (e.g. sample1_R1.fastq.gz and sample1_R2.fastq.gz)
#The script will ditch reads which have the overrepresented kmers present and their corresponding R1/R2 read (unlike the overrepresented_sequences script, the reverse comp of the kmer is not searched, because it is assumed they are pegged to the start or end of the read)
#This script assumes you have already ditched singletons using the overrepresented_sequences.R script.

library(stringi)
forwards <- list.files()[grep("R1",list.files())]
reverses <- list.files()[grep("R2",list.files())]
overrepresented <- as.matrix(read.table("overrepresented_kmers.fas"))

outputrecord <- matrix(NA, ncol=((dim(overrepresented)[1])/2)+3,nrow=length(forwards)+1)
outputrecord[1,] <- c("R1_name","R2_name","total_no_input_seq",overrepresented[seq(1,(dim(overrepresented)[1]),2),1])

for (j in seq(2,(dim(overrepresented)[1]),2)) {
overrepresented[j,1] <- paste("^",gsub("N","[[:alpha:]]",overrepresented[j,1],fixed=TRUE),sep="")
}

for (i in 1:length(forwards)) {
temp1 <- readLines(forwards[i])
temp2 <- readLines(reverses[i])
fastalines <- seq(2,length(temp1),4)
outputrecord[(i+1),1] <- forwards[i]
outputrecord[(i+1),2] <- reverses[i]
outputrecord[(i+1),3] <- length(fastalines)

for (j in seq(2,(dim(overrepresented)[1]),2)) {
tempcount <- which(((grepl(overrepresented[j,1],temp1[fastalines])) | (grepl(overrepresented[j,1],temp2[fastalines])))==TRUE)
outputrecord[(i+1),((j/2)+3)] <- length(tempcount)

if (length(tempcount)>0) {
tempcount <- sort(c(fastalines[tempcount],fastalines[tempcount]+1,fastalines[tempcount]+2,fastalines[tempcount]-1))
temp1 <- temp1[-tempcount]
temp2 <- temp2[-tempcount]
}
}

newFname <- gsub("R1","R1_kmerremoved",forwards[i])
newFname <- gsub(".gz","",newFname)
newRname <- gsub("R2","R2_kmerremoved",reverses[i])
newRname <- gsub(".gz","",newRname)

write.table(temp1,newFname,quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(temp2,newRname,quote=FALSE,col.names=FALSE,row.names=FALSE)

}

write.table(outputrecord,"outputrecord.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
