#!/bin/sh
# This script takes a bed-file of regions and uses heatmapper3.py to generate a signal matrix for each ChIP sample. Then, the signals are clustered according to distance and linkage methods indicated in the calling of the script.
# usage example: sh shellWrapper.sh	AllPeaks_merged		log2ratio	euclidean	ward	new
###########################################################################################################
#${1} bed-file with merged peak regions but _without_ .bed ending! e.g. RandSample_AllPeaks_merged - see below for hints on how to generate such a file
#${2} signal type, either one of: "log2ratio-noGC,"difference" (GC-corrected), "difference-InputScale","log2ratio-GC"
#${3} distance method; This must be one of "pearson", "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
#${4} linkage method; one of "ward", "single", "complete", "average", "mcquitty", "median" or "centroid"
#${5} whether to use an existing distance matrix (R-object) or create a new distance matrix: indicate Robj path or "new"
###########################################################################################################
## output files:
# peakMatrix_Complete_${1}_${2}.tab --> signal matrix
# peakMatrix_ChIPsignals_HMcollage_${1}_${2}.png --> visualizes individual signals for the regions taken into account for the signal matrix
# *dist.Robj --> R-object that contains the distances matrix
# *hclust.Robj --> R-object that contains the results of the hierarchical clustering
# dendrogram*png --> visualization of the hiearchical clustering
###########################################################################################################
## obtaining a bed-file of merged peaks - here's an example:
#cat ../../data/Tomek/peaks/MCRS1_merge_CAB.bed ../../data/Tomek/peaks/MOF_merge_CAB.bed ../../data/Tomek/peaks/MSL1_merge_CAB.bed ../../data/Tomek/peaks/MSL2_merge_CAB.bed ../../data/Tomek/peaks/NSL3_merge_CAB.bed > AllPeaks.bed
#/package/BEDTools-Version-2.12.0/bin/mergeBed -i AllPeaks.bed > AllPeaks_merged.bed
##########################################################################################################
### 1. signal matrix
for CHIP in MSL1 MSL2 MOF MCRS1 NSL3
do
		 /package/deepTools-1.5.7/bin/computeMatrix scale-regions  -S $SIGNAL -R ${1}.bed -F bigwig  -w 50 -m 1500 --regionsLabel "peak regions" --sortRegions no --colorMap Spectral_r --missingDataAsZero -T "${CHIP}merge_mlES unclustered ${2} signal" --startLabel "peak start" --endLabel "peak end" --outFileName peakMatrix_${1}_${2}_${CHIP}merge_mlES.png --outFileNameMatrix peakMatrix_${1}_${2}_${CHIP}merge_mlES.tab
done
#
montage peakMatrix_${1}_${2}_MOFmerge_mlES.png peakMatrix_${1}_${2}_MSL1merge_mlES.png peakMatrix_${1}_${2}_MSL2merge_mlES.png peakMatrix_${1}_${2}_MCRS1merge_mlES.png peakMatrix_${1}_${2}_NSL3merge_mlES.png -tile x1 -geometry +0.1+0.1 peakMatrix_ChIPsignals_HMcollage_${1}_${2}.png
paste peakMatrix_${1}_${2}_MOFmerge_mlES.tab peakMatrix_${1}_${2}_MSL1merge_mlES.tab peakMatrix_${1}_${2}_MSL2merge_mlES.tab peakMatrix_${1}_${2}_MCRS1merge_mlES.tab peakMatrix_${1}_${2}_NSL3merge_mlES.tab > peakMatrix_ALLChIPs_${1}_${2}.tab
echo "A signal matrix was produced: peakMatrix_ALLChIPs_${1}_${2}.tab.
The corresponding signal heatmaps are merged in peakMatrix_ChIPsignals_HMcollage.png"
rm peakMatrix_${1}_${2}_MOFmerge_mlES.tab peakMatrix_${1}_${2}_MSL1merge_mlES.tab peakMatrix_${1}_${2}_MSL2merge_mlES.tab peakMatrix_${1}_${2}_MCRS1merge_mlES.tab peakMatrix_${1}_${2}_NSL3merge_mlES.tab
rm peakMatrix_${1}_${2}_MOFmerge_mlES.png peakMatrix_${1}_${2}_MSL1merge_mlES.png peakMatrix_${1}_${2}_MSL2merge_mlES.png peakMatrix_${1}_${2}_MCRS1merge_mlES.png peakMatrix_${1}_${2}_NSL3merge_mlES.png
#
### 2. clustering
echo "Clustering of peakMatrix_ALLChIPs_${1}_${2}.tab with distance measure = ${3} and linkage method = ${4}"
/data/projects/muehlpfordt/scripts/Cluster1_computeHierarchicalClustering_FM_2013-03-11.R peakMatrix_ALLChIPs_${1}_${2}.tab ${3} ${4} ${2} ${5}
#
### 3. pruning dendrogram & heatmaps
echo "Dendrogram is pruned to 3, 4, 5 & 6 clusters"
for i in 3 4 5 6
do
/data/projects/muehlpfordt/scripts/Cluster2a_split_in_groups_plot_dendrogram_FM_outputUnsorted.R peakMatrix_ALLChIPs_${1}_${2}.tab hclust_${2}_${3}_${4}.Robj ${1}.bed ${i} clusters_${i}.bed dendrogram_${i}.pdf summary
/package/BEDTools/bin/bedtools shuffle -i ${1}.bed -g /data/misc/genomes/Mm/mm9/sizes/mm9_no_random.genome | head -2500 >> clusters_${i}.bed 

for CHIP in MOF MSL1 MSL2 MCRS1 NSL3
do
case ${2} in
		difference)
		SIGNAL=/data/projects/akhtar/bigwig/Tomek_mm9/differenceFiles_GCcorrected/difference_${CHIP}merge_mlES.bw;;
#
		log2ratio-noGC)
		SIGNAL=/data/projects/akhtar/bigwig/Tomek_mm9/log2ratio_${CHIP}merge_mlES_noGCcorrect.bw;;
#
	 	log2ratio-GC)
		SIGNAL=/data/projects/akhtar/bigwig/Tomek_mm9/log2ratio_GCcorrected/l2r_${CHIP}merge_mlES.bw;;
		
		difference-InputScale)
		SIGNAL=/data/projects/muehlpfordt/2013_Tomek/2013-03-08/difference_scaledInput/difference_${CHIP}merge_mlES_InputScaledUp.bw;;
esac
python /data/projects/ramirez/tools/heatmapper/heatmapper3.py reference-point -w 50 -a 2000 -b 2000 -S $SIGNAL -R clusters_${i}.bed -F bigwig --regionsLabel "random" --referencePoint center --refPointLabel "peak center" --colorMap RdYlBu --missingDataAsZero -T "${i} Clusters: ${CHIP}merge_mlES signal" --outFileName ${i}Clusters_merge_mlES_${CHIP}.png --sortRegions descend --whatToShow "heatmap and colorbar"
done
sh /data/projects/muehlpfordt/scripts/montage.sh ${i}Clusters_merge_mlES
done
echo "Heatmaps were produced for 3, 4, 5, 6 clusters" 
