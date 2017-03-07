library(stringr)

key <- as.matrix(read.table("key",sep="\t"))

fastafile <- as.matrix(read.table("name",sep="\t",quote=""))

outputsequence <- NULL

tempfile <- as.matrix(read.table(fastafile,sep="\t"))
ucelocus <- gsub(".fasta","",fastafile)

for (j in 1:dim(key)[1]) {
tempfile[(which(tempfile[,1]==key[j,1])),1] <- eval(parse(text=key[j,2]))
}

linestoditch <- matrix(FALSE,ncol=1,nrow=(dim(tempfile)[1]))

for (j in seq(2, dim(tempfile)[1],2)) {
linestoditch[j,1] <- nchar(gsub("-","",tempfile[j,1]))==0
if (linestoditch[j,1]==TRUE) {
linestoditch[(j-1),1] <- TRUE
}
}

removingmissing <- as.matrix(tempfile[(which(linestoditch[,1]==FALSE)),],ncol=1)

for (j in seq(2, dim(removingmissing)[1],2)) {
removingmissing[j,1] <- gsub("-","",removingmissing[j,1])
}

outputsequence <- rbind(outputsequence,removingmissing)

write.table(outputsequence,"temp",quote=FALSE,row.names=FALSE,col.names=FALSE)
