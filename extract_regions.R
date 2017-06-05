# Expect the feature table to be tab-delimited

working_dir <- "C:\\Users\\Alana\\Dropbox\\beetles\\mitogenomes_18S_28S\\mitogenome\\final"
fasta_suffix <- ".fasta"
feature_suffix <- "_feature_table.txt"

setwd(working_dir)

for (i in list.files(pattern=fasta_suffix)) {
    temp_fasta <-  read.table(i,header=FALSE,stringsAsFactors=FALSE,sep="\t")
    fasta <- gsub(":","",temp_fasta[1,1],fixed=TRUE)
    sequencepaste <- NULL
    for (j in 2:(dim(temp_fasta)[1])) {
      if ((length(grep(">",temp_fasta[j,1])))>0) {
      to_write <- toupper(sequencepaste)
      fasta <- rbind(fasta,to_write)
      fasta <- rbind(fasta,gsub(":","",temp_fasta[j,1],fixed=TRUE))
      sequencepaste <- NULL
      } else {
      sequencepaste <- paste(sequencepaste,temp_fasta[j,1],sep="")
      }
    }
    fasta <- rbind(fasta,toupper(sequencepaste))
    
    featurename <- gsub(fasta_suffix,feature_suffix,i,fixed=TRUE)
    feature <- readLines(featurename)
    for (j in 1:(length(feature))) {
      if ((length(grep(">",feature[j])))>0) {
          which_fasta <- which(fasta[,1]==gsub("Feature ","",feature[j]))
          fastaheader <- fasta[which_fasta,1]
          fastaseq <- unlist(strsplit(fasta[(which_fasta+1),1],""))
      } else {
          temp <- unlist(strsplit(feature[j],"\t"))
          if (!(is.na(suppressWarnings(as.numeric(temp[1])*as.numeric(temp[2]))))) {
              startpos <- as.numeric(temp[1])
              endpos <- as.numeric(temp[2])
          } else {
             if ((length(grep("gene",feature[j])))>0) {
                 genename <- unlist(strsplit(feature[j],"gene\t"))[2]
                 genename <- gsub("(","",genename,fixed=TRUE)
                 genename <- gsub(")","",genename,fixed=TRUE)
                 genename <- gsub("-.*","",genename,fixed=FALSE)
                 genename <- gsub("_.*","",genename,fixed=FALSE)
                 genename <- paste(genename,".fasta",sep="")
                 if (startpos < 0 ) {
                    geneseq <- paste(fastaseq[0:endpos],collapse="")
                    Ns <- paste(rep("N",(0-startpos)),collapse="")
                    geneseq <- paste(Ns,geneseq,sep="")                     
                 } else {  
                    geneseq <- paste(fastaseq[startpos:endpos],collapse="")
                 }    
                 write.table(fastaheader,genename,quote=FALSE, row.names=FALSE,col.names=FALSE,append=TRUE)
                 write.table(geneseq,genename,quote=FALSE, row.names=FALSE,col.names=FALSE,append=TRUE)
              }
           }
       }
    }
 }
