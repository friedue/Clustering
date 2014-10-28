### This script gives an overview of the different steps for the shell-
### based clustering routine using deepTools computeMatrix and the
### hclust function from R.
#######################
## Analysis-specific parameters
#######################
infolder=###!!! in this example, you could either indicate infolder=data or leave the infolder variable completely as it is clear from the folder structure where the input data can be found
outfolder=###!!! same as for infolder, here it might make sense to indicate outfolder=results, but then again, this is so trivial that it might be left out
scriptfolder=/data/projects/akhtar/bhardwaj/project1/Gitrepo/Clustering
rscript=/package/R-3.1.0/bin/Rscript ####!!!  this should be in your bashrc, not the script
bedtools=/package/BEDTools/bin/bedtools ####!!!  this should also be in your bashrc, not in the script
distanceMeasure=euclidean
linkage=ward
analysis= ###!!! this should be a unique identifier connected
## to the individual analysis, e.g. a combination of the linkage,
## distance and samples used or the like --> this should also be part
## of the output folder name (e.g. 2014-10_EuclidWard)
##########################
## The workflow
##########################
mkdir ${outfolder}

## compute individual data value matrices
for SAMPLE in MSL1_mlES MSL2_mlNPC
do
echo "computing the values for ${SAMPLE} and <REGIONFILE>"
sh ${scriptfolder}/a_generateMatrices.sh <SIGNALFILE> <REGIONFILE> ${outfolder} ${SAMPLE}
done

## rank transformation of each matrix individually,
## generating one composite matrix that will contain all matrices
## pasted next to each other

ls | grep "tab" > filenames.txt ####!!!! this needs to be more robust
echo "rank transforming each matrix"
${rscript} ${scriptfolder}/b_matrixProcessing.R \
filenames.txt \
${outfolder}/processedMatrix_${analysis}

## hierarchical clustering
echo "performing the hierachical clustering"
${rscript} ${scriptfolder}/c_HierarchicalClustering.R \
${outfolder}/processedMatrix_${analysis}.tab \ # input file for hclust
${distanceMeasure} ${linkage} new # various additional parameters for the R script

## prune the dendrogram at different levels to
## create clusters and produce bed files with clusters as output
for i in 2 3 4 5 6 7 8 9 10 
do
echo "dendrogram pruning for ${i} clusters"
${rscript} ${scriptfolder}/d_dendrogramPruning_outputUnsorted.R \
${outfolder}/processedMatrix_${analysis}.tab\
hclust_${analysis}.Robj \
 <REGIONFILE> \ ##!!!!!!!
${i} ${outfolder}/clusters_${i}.bed\
${outfolder}/dendrogram_${i}.pdf \
summary ###!!!! note that you can use this script with "individual"
##!!!! instead of "summary" to get the individual bed files per cluster!
# merge a shuffled region as control to the output bed files 
${bedtools} shuffle -i <REGIONFILE>\
-g /data/projects/misc/genomes/Mm/mm9/sizes/mm9_no_random.genome | \
head -2500 >> ${outfolder}/clusters_${i}.bed
echo "#random" >> ${outfolder}/linkagemethod_ward/clusters_${i}.bed
done


###!!!!! add the multiheatmapper version of computeMatrix and heatmapper
###!!!!! routine for the individual clusters


