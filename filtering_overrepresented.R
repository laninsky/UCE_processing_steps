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
blankboth <- sort(c(fastalines[blankboth],fastalines[blankboth]+1,fastalines[blankboth]+2,fastalines[blankboth]-1))
temp1 <- temp1[-blankboth]
temp2 <- temp2[-blankboth]


blankreverses <- which((temp2[fastalines]=="")==TRUE)




temp1[fastalines[blankforwards]]

which(grepl(temp1[fastalines],"")==TRUE)

singletonout <- 
