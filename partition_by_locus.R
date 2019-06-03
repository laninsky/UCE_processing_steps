#1. Getting a list of the phylip files in the phylip_w_missing folder
listofphylipfiles <- list.files(pattern="*.phy*")

#2. Reading in the first file and using this to set up the variables 
# we'll use to record aspects of the data
locusname <- gsub("-","_",gsub(".phy.*","",listofphylipfiles[1]))
temp <- readLines(listofphylipfiles[1])

notaxa <- as.numeric(unlist(strsplit(temp[1],"\\s"))[2])
finalnameseq <- matrix(ncol=2,nrow=notaxa)

firstrows <- unlist(strsplit(temp[2:(1+notaxa)],"\\s"))
if (length(which(firstrows==""))>0) {
  firstrows <- firstrows[-which(firstrows=="")]
}
finalnameseq[,1] <- unlist(strsplit(firstrows,"\\s"))[seq(1,length(firstrows),round(length(firstrows)/notaxa))]
restofseq <- gsub("\\s","",unlist(strsplit(firstrows,"\\s"))[-seq(1,length(firstrows),round(length(firstrows)/notaxa))])

for (i in 1:dim(finalnameseq)[1]) {
  finalnameseq[i,2] <- paste(restofseq[
    (1+(i-1)*length(restofseq)/notaxa):((1+(i-1)*length(restofseq))/notaxa+length(restofseq)/notaxa)],collapse="")
}

if (length(which(temp==""))>0) {
  temp <- temp[-which(temp=="")]
}

temp <- temp[-(1:(notaxa+1))]
temp <- gsub("\\s","",temp)
names(temp) <- rep(finalnameseq[,1],length(temp)/notaxa)

temptempseq <- matrix(nrow=dim(finalnameseq)[1])

for (i in 1:dim(finalnameseq)[1]) {
  temptempseq[i,1] <- paste(temp[which(names(temp)==finalnameseq[i,1])],collapse="")
  finalnameseq[i,2] <- paste(finalnameseq[i,2],temptempseq[i,1],collapse="")
}

finalnameseq <- gsub("\\s","",finalnameseq)

chrsettable <- paste("CHARSET ",locusname, " = ",1,"-",nchar(finalnameseq[1,2]),";",sep="")

#3. Pulling each of the remaining phylip files through
for (j in listofphylipfiles[-1]) {
  # Reading the uce locus
  temp <- readLines(j)
  # Obtaining the name
  locusname <- gsub("-","_",gsub(".phy.*","",j))
  # Extracting the sequence
  firstrows <- unlist(strsplit(temp[2:(1+notaxa)],"\\s"))
  if (length(which(firstrows==""))>0) {
    firstrows <- firstrows[-which(firstrows=="")]
  }
  
  tempseq <- matrix(nrow=dim(finalnameseq)[1],ncol=2)
  tempseq[,1] <- unlist(strsplit(firstrows,"\\s"))[seq(1,length(firstrows),round(length(firstrows)/notaxa))]
  restofseq <- gsub("\\s","",unlist(strsplit(firstrows,"\\s"))[-seq(1,length(firstrows),round(length(firstrows)/notaxa))])
  # Getting the start position for this locus
  starpos <- nchar(finalnameseq[1,2])+1
  # Getting the first rows of sequence together
  for (i in 1:dim(tempseq)[1]) {
    tempseq[i,2] <- paste(restofseq[
      (1+(i-1)*length(restofseq)/notaxa):((1+(i-1)*length(restofseq))/notaxa+length(restofseq)/notaxa)],collapse="")
  }
  if (length(which(temp==""))>0) {
    temp <- temp[-which(temp=="")]
  }
  temp <- temp[-(1:(notaxa+1))]
  temp <- gsub("\\s","",temp)
  names(temp) <- rep(tempseq[,1],round(length(temp)/notaxa))

  temptempseq <- matrix(nrow=dim(finalnameseq)[1])
  
  for (i in 1:dim(tempseq)[1]) {
    temptempseq[i,1] <- paste(temp[which(names(temp)==tempseq[i,1])],collapse="")
    tempseq[i,2] <- paste(tempseq[i,2],temptempseq[i,1],collapse="")
  }
  
  tempseq <- gsub("\\s","",tempseq)

  for (i in 1:dim(finalnameseq)[1]) {
    finalnameseq[i,2] <- paste(finalnameseq[i,2],tempseq[which(tempseq[,1]==finalnameseq[i,1]),2],collapse="")
  }

  finalnameseq <- gsub("\\s","",finalnameseq)
  
  tempchrsettable <- paste("CHARSET ",locusname, " = ",starpos,"-",nchar(finalnameseq[1,2]),";",sep="")
  chrsettable <- rbind(chrsettable,tempchrsettable)
  
}

#4. Writing out aspects of the nexus file
write("#NEXUS","outputdata.nex")
write("","outputdata.nex",append=TRUE)
write("BEGIN DATA;","outputdata.nex",append=TRUE)
write(paste("DIMENSIONS NTAX=",notaxa," NCHAR=",nchar(finalnameseq[1,2]),";",sep=""),"outputdata.nex",append=TRUE)
write("FORMAT DATATYPE=DNA GAP=- MISSING=?;","outputdata.nex",append=TRUE)
write("MATRIX","outputdata.nex",append=TRUE)
write("","outputdata.nex",append=TRUE)
write.table(finalnameseq,"outputdata.nex",append=TRUE,quote = FALSE,row.names = FALSE,col.names = FALSE)
write(";","outputdata.nex",append=TRUE)
write("","outputdata.nex",append=TRUE)
write("END;","outputdata.nex",append=TRUE)
write("","outputdata.nex",append=TRUE)
write("BEGIN SETS;","outputdata.nex",append=TRUE)
write("","outputdata.nex",append=TRUE)
write("[loci]","outputdata.nex",append=TRUE)
write.table(chrsettable,"outputdata.nex",append=TRUE,quote = FALSE,row.names = FALSE,col.names = FALSE)
write("","outputdata.nex",append=TRUE)

charpartition <- "CHARPARTITION loci = "
for (i in 1:dim(chrsettable)[1]) {
  temp <- gsub(" = .*","",gsub("CHARSET ","",chrsettable[i,1]))
  charpartition <- paste(charpartition,i,":",temp,", ",sep="")
}

charpartition <- paste(substr(charpartition,1,(nchar(charpartition)-2)),";",sep="")

write(charpartition,"outputdata.nex",append=TRUE)
write("","outputdata.nex",append=TRUE)
write("END;","outputdata.nex",append=TRUE)

