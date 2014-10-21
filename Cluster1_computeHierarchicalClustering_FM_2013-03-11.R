#! /usr/bin/Rscript --vanilla
# usage: ./computeHierarchicalClustering.R peakMatrix_ALL.tab euclidean ward log2ratio new
# args: [1] signal matrix for all peaks and all ChIPs
#       [2] distance method; This must be one of "pearson", "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
#       [3] linkage method; one of "ward", "single", "complete", "average", "mcquitty", "median" or "centroid"
#       [4] what kind of signal was used for signal matrix, e.g. log2ratio or difference
#       [5] either "new" or  an existing distance.Robj
library(GMD)
args <- commandArgs(TRUE)

### defining filenames for saving
outFile_dist = paste("dist_",args[4],"_", args[2], ".Robj", sep="")
outFile_hclust = paste("hclust_",args[4],"_",  args[2],"_", args[3], ".Robj", sep="")


if (args[5] != "new"){
  load(args[5])
  } else {
    ma = read.table(args[1])
    if (args[2] == "pearson"){
      #corPearson = cor(t(as.matrix(ma)), method="pearson", use="complete.obs")
      #distanceMM = as.dist((1 - corPearson)/2)
      #hh.dist = hclust(distanceMM, method = args[3])
      distanceMM = as.dist(1-cor(t(as.matrix(ma)), method='pearson') )
      } else {
        distanceMM=gdist(as.matrix(ma), method = args[2]) # 3013-03-10: changed dist() to gdist()
      }
    save(distanceMM, file=outFile_dist)
  }

if( (args[3] == "ward") || (args[3]=="median") || (args[3] == "centroid") ){
  hh.dist = hclust(distanceMM^2, method=args[3])
}else {
  hh.dist = hclust(distanceMM, method=args[3])
}
save(hh.dist, file = outFile_hclust)


# plotting a dendrogram of the hierarchical clustering
png(file=paste("dendrogram",args[4],args[2], args[3],"png", sep="."), width=1200, height=400)
plot(hh.dist, labels=FALSE, main=paste(args[4],args[2], args[3]))
dev.off()
