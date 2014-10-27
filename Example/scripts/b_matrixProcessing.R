#! /usr/bin/Rscript --vanilla
# usage: Rscript ListOfChIPMatrices Outname
# e.g. Rscript filelist.tab peakMatrix_processed

args <- commandArgs(TRUE)

FileList = as.vector(read.table(args[1])[,1])

### function for data processing
ProcessTabFileForCluster <- function(InFile, Method = "rank"){
  # Method can bei either "linear" or "rank"
  mm = as.matrix(read.table(InFile))
  if(Method == "linear"){
    mm_filt = mm/max(mm)
    mm_filt[mm_filt < 0] = 0
  }else{
    mm[mm < 0] = 0
    mm_filt = matrix(rank(mm), ncol = dim(mm)[2])
  }
  return(mm_filt)
}

### read in, process data, generate big list that stores all the processed data matrices
MatrixList = list()
for(FILE in FileList){
  MatrixList[[FILE]] = ProcessTabFileForCluster(InFile=FILE) # rank-process all ChIP matrices
}

### save processed data
output = data.frame(do.call(cbind, MatrixList)) # paste all ChIP matrices together
write.table(output,file =paste(args[2],".tab", sep = ""), col.names =F, row.names = F, quote = F, sep='\t')
