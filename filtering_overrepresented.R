#To run this script modify the path to "overrepresented" to your fasta file containing the overrepresented sequences across your samples
#The script should be placed in a folder with your F and R reads and executed by `Rscript filtering_overrepresentatives.R`
#In order for the script to work, you should have a R1 and R2 file for each individual (e.g. sample1_R1.fastq.gz and sample1_R2.fastq.gz)
#The script will ditch reads which have the overrepresented sequences present (or the reverse complement) as well as their corresponding R1/R2 read
#It will also ditch blank sequences and write the singleton out into a file prefixed with singleton

library(stringi)
forwards <- list.files()[grep("R1",list.files())]
reverses <- list.files()[grep("R2",list.files())]
overrepresented <- as.matrix(read.table("/home/a499a400/manu/combined/overrepresented_sequences.fas"))

outputrecord <- matrix(NA, ncol=((dim(overrepresented)[1])/2)+6,nrow=length(forwards)+1)
outputrecord[1,] <- c("R1_name","R2_name","total_no_input_seq","both_blank_seqs_removed","R1_only_removed","R2_only_removed",overrepresented[seq(1,(dim(overrepresented)[1]),2),1])

for (i in 1:length(forwards)) {
temp1 <- readLines(forwards[i])
temp2 <- readLines(reverses[i])
fastalines <- seq(2,length(temp1),4)
outputrecord[(i+1),1] <- forwards[i]
outputrecord[(i+1),2] <- reverses[i]
outputrecord[(i+1),3] <- length(fastalines)

blankboth <- which(((temp1[fastalines]=="") & (temp2[fastalines]==""))==TRUE)
outputrecord[(i+1),4] <- length(blankboth)

if (length(blankboth)>0) {
blankboth <- sort(c(fastalines[blankboth],fastalines[blankboth]+1,fastalines[blankboth]+2,fastalines[blankboth]-1))
temp1 <- temp1[-blankboth]
temp2 <- temp2[-blankboth]
}

for (j in seq(2,(dim(overrepresented)[1]),2)) {
revcomp <- stri_reverse(chartr("acgtACGT","tgcaTGCA",overrepresented[j,1]))

tempcount <- which(((grepl(overrepresented[j,1],temp1[fastalines])) | (grepl(overrepresented[j,1],temp2[fastalines])) | (grepl(revcomp,temp1[fastalines])) | (grepl(revcomp,temp2[fastalines])))==TRUE)

outputrecord[(i+1),((j/2)+6)] <- length(tempcount)

if (length(tempcount)>0) {
tempcount <- sort(c(fastalines[tempcount],fastalines[tempcount]+1,fastalines[tempcount]+2,fastalines[tempcount]-1))

temp1 <- temp1[-tempcount]
temp2 <- temp2[-tempcount]
}
}

blankforward <- which((temp1[fastalines]=="")==TRUE)
outputrecord[(i+1),5] <- length(blankforward)

if (length(blankforward)>0) {
blankforward <- sort(c(fastalines[blankforward],fastalines[blankforward]+1,fastalines[blankforward]+2,fastalines[blankforward]-1))
newname <- gsub("R1","singleton",forwards[i])
newname <- gsub(".gz","",newname)
write.table(temp2[blankforward],newname,quote=FALSE,col.names=FALSE,row.names=FALSE)

temp1 <- temp1[-blankforward]
temp2 <- temp2[-blankforward]
}

blankreverse <- which((temp2[fastalines]=="")==TRUE)
outputrecord[(i+1),6] <- length(blankreverse)

if (length(blankreverse)>0) {
blankreverse <- sort(c(fastalines[blankreverse],fastalines[blankreverse]+1,fastalines[blankreverse]+2,fastalines[blankreverse]-1))
newname <- gsub("R1","singleton",forwards[i])
newname <- gsub(".gz","",newname)
write.table(temp1[blankreverse],newname,quote=FALSE,append=TRUE,col.names=FALSE,row.names=FALSE)

temp1 <- temp1[-blankreverse]
temp2 <- temp2[-blankreverse]
}

newFname <- gsub("R1","R1_repeatremoved",forwards[i])
newFname <- gsub(".gz","",newFname)
newRname <- gsub("R2","R2_repeatremoved",reverses[i])
newRname <- gsub(".gz","",newRname)

write.table(temp1,newFname,quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(temp2,newRname,quote=FALSE,col.names=FALSE,row.names=FALSE)

}

write.table(outputrecord,"outputrecord.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
