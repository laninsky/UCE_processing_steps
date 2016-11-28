k <- 1
partitionfile <- NULL
seqblock <- NULL
for (i in list.files(pattern="concatphylip")) {
temp <- as.matrix(read.table(paste(i,"/",(list.files(path=i)),sep="")))
temppartition <- paste("DNA, ",i,"=",k,"-",(k+as.numeric(temp[1,2])-1),sep="")
partitionfile <- rbind(partitionfile,temppartition)
k <- k+as.numeric(temp[1,2])
temp <- temp[-1,]
temp <- temp[order(temp[,1]),]
seqblock <- as.matrix(paste(seqblock,temp[,2],sep=""),ncol=1)
}
temp <- as.matrix(read.table(paste(i,"/",(list.files(path=i)),sep="")))
notaxa <- temp[1,1]
temp <- temp[-1,]
temp <- temp[order(temp[,1]),]
datamatrix <- cbind(temp[,1],seqblock)
titleline <- paste(notaxa,(k-1))

finaldatamatrix <- matrix(NA,ncol=1,nrow=as.numeric(notaxa))
for (i in 1:dim(datamatrix)[1]) {
finaldatamatrix[i,1] <- paste(datamatrix[i,1],paste(rep(" ",(10-nchar(datamatrix[i,1]))),collapse=""),datamatrix[i,2],sep="")
}

finaldatamatrix <- rbind(titleline,finaldatamatrix)

write.table(finaldatamatrix,"concatphylip.phylip",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(partitionfile,"partitions.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
