data <- as.matrix(read.table("C://Users//a499a400//Dropbox//beetles//stephen//50perc_w_missing_no_title.txt"))

data <- data[-1,]

output <- NULL

for (i in 1:(dim(data)[1])) {
    output <- rbind(output,(unlist(strsplit(data[i,2],""))))
}

missing_by_sample <- c("samples","mean")
for (i in 1:(dim(data)[1])) {
   tempmissing <- c(data[i,1],((sum(output[i,1:(dim(output)[2])]=="?")/(dim(output)[2]))*100))
   missing_by_sample <- rbind(missing_by_sample,tempmissing)
}

nomissing <- NULL

for (i in 1:(dim(output)[2])) {
    nomissing <- c(nomissing,sum(output[1:(dim(output)[1]),i]=="?"))
}

output <- c("number of missing samples","number of sites")
for (i in 0:(dim(data)[1])) {
    tempoutput <- c(i,sum(nomissing==i))
    output <- rbind(output,tempoutput)
}
