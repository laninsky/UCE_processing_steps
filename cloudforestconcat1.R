models <- as.matrix(read.table("output_models.txt",sep="\t"))
models[,2] <- gsub("AICc-","",models[,2])
modelnames <- unique(models[,2])

modelnamecommands <- paste("mkdir ",modelnames,sep="")
modelnamecommands <- as.matrix(modelnamecommands,ncol=1,nrow=(length(modelnames)))

copycommands <- paste("mv ",models[,1],".phylip ", models[,2],sep="")
copycommands <- as.matrix(copycommands,ncol=1,nrow=(dim(models)[1]))

phylucecommands <- paste("phyluce_align_convert_one_align_to_another --alignments ",modelnames," --output ",modelnames,".nexus --input-format phylip --output-format nexus --cores 4; phyluce_align_format_nexus_files_for_raxml --alignments ",modelnames,".nexus --output ",modelnames,".concatphylip",sep="")
phylucecommands <- as.matrix(phylucecommands,ncol=1,nrow=(length(modelnames)))

commands <- rbind(modelnamecommands,copycommands,phylucecommands)

write.table(commands,"concat_by_model.sh",quote=FALSE,row.names=FALSE,col.names=FALSE)
