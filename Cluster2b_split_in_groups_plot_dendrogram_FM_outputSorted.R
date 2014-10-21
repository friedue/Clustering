#! /usr/bin/Rscript --vanilla
# example usage: ./split_in_groups_plot_dendrogram.R peakMatrix_ALL.tab hclust_Euclid.Robj AllPeaks_merged.bed 3 clusters_3.bed dendrogram_3.pdf
# args: [1] input: matrix with all signals for all peaks
#       [2] input: hclust result
#       [3] input: merged peaks bed-file
#       [4] number of clusters
#       [5] output name .bed file
#       [6] output name dendrograms
#       [7] indicate "individual" if individual bed-files should be generated for each cluster (i.e. will result in 3 bed-files when 3 clustes are chosen)
library(ggdendro)
library(maptree)

args <- commandArgs(TRUE)

ma = read.table(args[1])
load(args[2]) # this loads an object called hh.dist = distance matrix with hierarchical clustering

cc = clip.clust(hh.dist, data=ma, k=as.integer( args[4] )) # Reduces a hierarchical cluster tree to a smaller tree either by pruning until a given number of observation groups remain, or by pruning tree splits below a given height

### saving information which peak region belongs to which cluster
bed = read.table(args[3]) # merged peaks bed-file
outFile=args[5] 
LENGTHS=as.numeric()
write.table(file="hclust_ordered.bed", bed[hh.dist$order,],sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
for(i in 1:length(cc$membership)) {
  print(paste("size of cluster",i, "is", length(cc$membership[[i]])))
  LENGTHS[i]=length(cc$membership[[i]])
  
  if (args[7]=="individual"){
    write.table(file=paste(outFile,i, sep="_"), bed[cc$membership[[i]],],sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE) # for saving individual bed-files per cluster-group
    } else {
      if (i<2){
        write.table(file=outFile, append=TRUE,bed[hh.dist$order[1:(length(cc$membership[[i]])-1)],],sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
        write.table(file=outFile, paste("#Cluster_", i, sep=''), append=TRUE, sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
      }else{
        j=i-1
        START=sum(LENGTHS[1:j]) + 1
        END=sum(LENGTHS[1:i]) - 1
        write.table(file=outFile, append=TRUE,bed[hh.dist$order[START:END],],sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
        write.table(file=outFile, paste("#Cluster_", i, sep=''), append=TRUE, sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE) # inserts hashed-comment into bed-file so heatmapper will recognize the different groups
      }
    }
}

### generating dendrograms
if(as.integer(args[4])>2) {
	tiff(file=paste("with_Totals_", args[6],sep=""), height= 300, width=500, type="cairo")
	draw.clust(cc, nodeinfo=TRUE)
	dev.off()
}
dev.off()
