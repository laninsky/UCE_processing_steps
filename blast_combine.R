#Reading in the contigs > 1000bp
contig_file_name <- as.matrix(read.table("contig_file_name"))
contigs <- as.matrix(read.table(contig_file_name[1,1]))

#Reading in the blast results
blast <- as.matrix(read.table("blast_output.txt"))

#Lines in the blast table where the qstart or qend is within 50bp of the beginning of the query
qbeginning <- unique(c(which(as.numeric(blast[,6])<=50)),(which(as.numeric(blast[,7])<=50)))
#Lines in the blast table where the qstart or qend is within 50bp of the end of the query
qend <- unique(c(which(as.numeric(blast[,6])>=(as.numeric(blast[,5])-50))),(which(as.numeric(blast[,7])>=(as.numeric(blast[,5])-50))))
#Lines in the blast table where the sstart or send is within 50bp of the begining of the subject
sbeginning <- unique(c(which(as.numeric(blast[,9])<=50)),(which(as.numeric(blast[,10])<=50)))
#Lines in the blast table where the sstart or send is within 50bp of the end of the subject
send <- unique(c(which(as.numeric(blast[,9])>=(as.numeric(blast[,8])-50))),(which(as.numeric(blast[,10])>=(as.numeric(blast[,8])-50))))

#Lines in the blast table where the blast match is at the beginning of the query, and the end of the subject
qbeg_sbeg <- qbeginning[qbeginning %in% sbeginning]
#Lines in the blast table where the blast match is at the beginning of the query, and the beginning of the subject
qbeg_send <- qbeginning[qbeginning %in% send]
#Lines in the blast table where the blast match is at the end of the query, and the beginning of the subject
qend_sbeg <- qbeginning[qend %in% sbeginning]
#Lines in the blast table where the blast match is at the end of the query, and the end of the subject
qend_send <- qbeginning[qend %in% send]

#Subsetting these promising rows (overhangs won't be in the middle of contigs)
subsetrows <- unique(c(qbeg_sbeg,qbeg_send,qend_sbeg,qend_send))
blast <- blast[subsetrows,]

#1     2     3       4     5     6     7     8     9    10
#qacc sacc pident length qlen qstart qend slen sstart send evalue bitscore

# Stepping through the blast table
for(i in 1:(dim(blast)[1])) {
#If either sequence name is na, skipping to the next row
  if(is.na(blast[i,1]) | is.na(blast[i,2])) {
    next
  }
  #Pulling out the sequence from the contig file corresponding to the blast table
  #Recording the row numbers just in case we need to delete these rows later on (if we merge sequence)
  temp1_name_row <- which(contigs[,1] %in% paste(">",blast[i,1],sep=""))
  temp2_name_row <- which(contigs[,1] %in% paste(">",blast[i,2],sep=""))
  temp1_seq_row <- temp1_name_row+1
  temp2_seq_row <- temp2_name_row+1
  temp1_seq <- contigs[temp1_seq_row,1]
  temp2_seq <- contigs[temp2_seq_row,1]
  
  #Figuring out where the match is on the two contigs
  if( (as.numeric(blast[i,6])<=50) | (as.numeric(blast[i,7])<=50) ) {
    temp1_stat <- "begin"
    } else {
    temp1_stat <- "end"
  }
  if( (as.numeric(blast[i,9])<=50) | (as.numeric(blast[i,10])<=50) ) {
    temp2_stat <- "begin"
    } else {
    temp2_stat <- "end"
  }
  if(temp1_stat=="begin") {
    if(temp2_stat=="begin") {
    #If both matches are at the beginning of the contigs, one of the contigs should be reverse complemented if it is an overhang
    #If not, we skip to the next line of the blast table
      if((((as.numeric(blast[i,7])-as.numeric(blast[i,6]))*(as.numeric(blast[i,10])-as.numeric(blast[i,9]))))>0) {
        next
      }
      #Otherwise we are going to reverse comp the first sequence in the table, after removing the overhanging bases
      outputseq <- unlist(strsplit(temp1_seq,""))
      if(as.numeric(blast[i,6])>as.numeric(blast[i,7])) {
        cut <- as.numeric(blast[i,6])
      } else {
        cut <- as.numeric(blast[i,7])
      }
      outputseq <- outputseq[-1:-cut]
      outputseq <- paste(outputseq,collapse="")
      outputseq <- toString(reverseComplement(DNAString(outputseq)))
      # We then paste the reversecomped sequence to the second
      outputseq <- paste(outputseq,temp2_seq,collapse="")
    } else {
    #This section is for matches at the beginning of contig 1, and at the end of contig 2. Neither sequence should be reverse comped in this case.
      if((((as.numeric(blast[i,7])-as.numeric(blast[i,6]))*(as.numeric(blast[i,10])-as.numeric(blast[i,9]))))<0) {
        next
      }
      #Here we are snipping off the overlapping region from contig 1 so we can stick them together
      outputseq <- unlist(strsplit(temp1_seq,""))
      if(as.numeric(blast[i,6])>as.numeric(blast[i,7])) {
        cut <- as.numeric(blast[i,6])
      } else {
        cut <- as.numeric(blast[i,7])
      }
      outputseq <- outputseq[-1:-cut]
      outputseq <- paste(outputseq,collapse="")
      #We stick them together with contig 2 first, then contig 1
      outputseq <- paste(temp2_seq,outputseq,collapse="")
  } else {
    # 
    if(temp2_stat=="begin") {
      if((((as.numeric(blast[i,7])-as.numeric(blast[i,6]))*(as.numeric(blast[i,10])-as.numeric(blast[i,9]))))<0) {
        next
      }
      outputseq <- unlist(strsplit(temp2_seq,""))
      if(as.numeric(blast[i,10])>as.numeric(blast[i,9])) {
        cut <- as.numeric(blast[i,10])
      } else {
        cut <- as.numeric(blast[i,9])
      }
      outputseq <- outputseq[-1:-cut]
      outputseq <- paste(outputseq,collapse="")
      outputseq <- paste(temp1_seq,outputseq,collapse="")
    } else {      
        if((((as.numeric(blast[i,7])-as.numeric(blast[i,6]))*(as.numeric(blast[i,10])-as.numeric(blast[i,9]))))>0) {
          next
        }
  
  
  # code to remove these sequences from the rest of the table #
  rm_rows <- c(temp1_name_row,temp2_name_row,temp1_seq_row,temp2_seq_row)
  contigs <- contigs[-rm_rows,]
  
  

title_lines <- which(contigs[,1] %in% contigs_names)

output <- NULL
for(i in title_lines) {
  output <- rbind(output,contigs[i,1])
  output <- rbind(output,contigs[(i+1),1])
  for(j in (i+2):(dim(contigs)[1])) {
    if(grepl(">",contigs[j,1],fixed=TRUE)) {
      break
    }
    output[(dim(output)[1]),1] <- paste(output[(dim(output)[1]),1],contigs[j,1],sep="")
  }
}

names <- unlist(strsplit(contig_file_name,"/"))
names <- unlist(strsplit(names[length(names)],"_"))
names <- paste(names[1],".fasta",sep="")

write.table(output,names,quote=FALSE,row.names=FALSE,col.names=FALSE)
