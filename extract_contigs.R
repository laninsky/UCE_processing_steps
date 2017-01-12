contigs_names <- as.matrix(read.table("extract_contigs.txt"))
contig_file_name <- as.matrix(read.table("contig_file_name"))
contigs <- as.matrix(read.table(contig_file_name[1,1]))

contigs_names <- paste(">",contigs_names,sep="")

title_lines <- which(contigs[,1] %in% contigs_names)

output <- NULL
for(i in title_lines) {
  output <- rbind(output,contigs[i,1])
  for(j in (i+1):(dim(contigs)[1])) {
    if(grepl(">",contigs[j,1],fixed=TRUE)) {
      break
    }
    output <- rbind(output,contigs[j,1])
  }
}

names <- unlist(strsplit(contig_file_name,"/"))
names <- unlist(strsplit(names[length(names)],"_"))
names <- paste(names[1],".fasta",sep="")

write.table(output,names,quote=FALSE,row.names=FALSE,col.names=FALSE)
