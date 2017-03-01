missing_data <- function(data_file) {

data <- as.matrix(read.table(data_file))

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

write.table(missing_by_sample,"missing_by_sample.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)    
    
nomissing <- NULL

for (i in 1:(dim(output)[2])) {
    nomissing <- c(nomissing,sum(output[1:(dim(output)[1]),i]=="?"))
}
    
output <- c("number of missing samples","number of sites")
for (i in 0:(dim(data)[1])) {
    tempoutput <- c(i,sum(nomissing==i))
    output <- rbind(output,tempoutput)
}
    
write.table(output,"missing_by_site.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
    
}
