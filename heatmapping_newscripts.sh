
## Vivek, 21-22rd October 2014.. 
## The modified scripts for the chIP results clustering for heatmap generation
## Parameters for Clustering: distance=euclidean, linkage method= ward, existing matrix= "new", signaltype= log2ratio
# 

infolder=/data/projects/akhtar/shared/project1/Tomek_Tasneem_Matthew_2014_AlleleSpecific
outfolder=/data/projects/akhtar/bhardwaj/project1
scriptfolder=/data/projects/akhtar/bhardwaj/project1/Gitrepo/Clustering
rscript=/package/R-3.1.0/bin/Rscript
bedtools=/package/BEDTools/bin/bedtools

## 1.Running computeMatrix on one file 
# A) fem ES cells

for CHIP in MSL1 MSL2 MOF MCRS1 NSL3; 
do CELL=femES; 
/package/deepTools-1.5.9.1/bin/computeMatrix scale-regions -S ${infolder}/bigWigs/${CHIP}merge-Inputmerge_${CELL}_l2r_readCount_optimalAlign.bw \
-R ${infolder}/Peaks/bedfiles/Peaks_ES-NPC_mergedByBedtools.bed \
-bs 50 -m 1500 --sortRegions no --missingDataAsZero \
--outFileNameMatrix ${outfolder}/peakMatrix_${CHIP}_${CELL}.tab --outFileName peakMatrix_outfile_${CHIP}_${CELL};
done

# B) femNPC's

for CHIP in MSL2 MOF NSL3; 
do CELL=femNPC; 
/package/deepTools-1.5.9.1/bin/computeMatrix scale-regions -S ${infolder}/bigWigs/${CHIP}merge-Inputmerge_${CELL}_l2r_readCount_optimalAlign.bw -R ${infolder}/Peaks/bedfiles/Peaks_ES-NPC_mergedByBedtools.bed -bs 50 -m 1500 --sortRegions no --missingDataAsZero --outFileNameMatrix ${outfolder}/peakMatrix_${CHIP}_${CELL}.tab --outFileName ${outfolder}/peakMatrix_outfile_${CHIP}_${CELL}; 
done

# The result is a matrix with signal values for the heatmapping (saved in /bhardwaj/projects1 folder)


## 2. Running Clustering with pre-prepared R scripts

## A) Signal matrix rank transformation (in /bhardwaj/projects1 folder)

ls | grep "tab" > filenames.txt
${rscript} ${scriptfolder}/Cluster0_matrixProcessing.R filenames.txt Cluster0_output_ES-NPC

## B) Clustering (hclust): parameters for this run: distance=euclidean, linkage method= complete, existing matrix= "new", signaltype= log2ratio
# for linkage method= complete, change the name of outfolder to ${outfolder}/linkagemethod_complete
linkage=complete	
${rscript} ${scriptfolder}/Cluster1_computeHierarchicalClustering_FM_2013-03-11.R ${outfolder}/Cluster0_output_ES-NPC.tab euclidean ${linkage} log2ratio-noGC new	

# Create clusters and produce bed files with clusters as output
for i in 7 8 9 10 
do
${rscript} ${scriptfolder}/Cluster2a_split_in_groups_plot_dendrogram_FM_outputUnsorted.R ${outfolder}/peakMatrix_ALLChIPs_ES-NPC.tab hclust_log2ratio-noGC_euclidean_${linkage}.Robj ${outfolder}/Peaks_ES-NPC_mergedByBedtools.bed ${i} ${outfolder}/linkagemethod_ward/clusters_${i}.bed ${outfolder}/linkagemethod_ward/dendrogram_${i}.pdf summary
#merge a shuffled region as control to the output bed files 
${bedtools} shuffle -i ${outfolder}/Peaks_ES-NPC_mergedByBedtools.bed -g /data/projects/misc/genomes/Mm/mm9/sizes/mm9_no_random.genome | 
head -2500 >> ${outfolder}/linkagemethod_ward/clusters_${i}.bed
echo "#random" >> ${outfolder}/linkagemethod_ward/clusters_${i}.bed
done


## 3) ComputeMatrix with the new clustered bed files 
#A) For femES cells 
for i in 5
do
for CHIP in MOF MSL1 MSL2 MCRS1 NSL3
do
echo "preparing matrix for ${CHIP} with ${i} clusters"
/package/deepTools-1.5.9.1/bin/computeMatrix reference-point -bs 50 -a 2000 -b 2000 \
-S ${infolder}/bigWigs/${CHIP}merge-Inputmerge_femES_l2r_readCount_optimalAlign.bw \
-R ${outfolder}/linkagemethod_complete/clusters_${i}.bed \
--referencePoint center --missingDataAsZero --skipZeros \
--outFileName ${outfolder}/linkagemethod_complete/peakMatrix_${i}_clustered_outfile_${CHIP}_femES.tab.gz;
done
done

#Then finally generate heatmaps using heatmapper and merge them togather 
i=5
for i in 5
do
for CHIP in MOF MSL1 MSL2 MCRS1 NSL3
do
echo "preparing heatmap for ${CHIP} with ${i} clusters"
/package/deepTools-1.5.9.1/bin/heatmapper -m ${outfolder}/linkagemethod_complete/peakMatrix_${i}_clustered_outfile_${CHIP}_femES.tab.gz --refPointLabel "peak center" --colorMap PuBuGn -T "${i} Clusters: ${CHIP}merge_femES signal" --outFileName ${outfolder}/linkagemethod_complete/heatmaps/${i}Clusters_merge_femES_${CHIP}.png --whatToShow "heatmap and colorbar" 	
done
done

# B) For femNPC cells:

for i in 5
do
for CHIP in MSL2 MOF NSL3
do
echo "preparing matrix for ${CHIP} with ${i} clusters"
/package/deepTools-1.5.9.1/bin/computeMatrix reference-point -bs 50 -a 2000 -b 2000 \
-S ${infolder}/bigWigs/${CHIP}merge-Inputmerge_femNPC_l2r_readCount_optimalAlign.bw \
-R ${outfolder}/linkagemethod_complete/clusters_${i}.bed \
--referencePoint center --missingDataAsZero --skipZeros \
--outFileName ${outfolder}/linkagemethod_complete/peakMatrix_${i}_clustered_outfile_${CHIP}_femNPC.tab.gz;
done
done

#heatmap
for i in 5
do
for CHIP in MOF MSL2 NSL3
do
echo "preparing heatmap for ${CHIP} with ${i} clusters"
/package/deepTools-1.5.9.1/bin/heatmapper -m ${outfolder}/linkagemethod_complete/peakMatrix_${i}_clustered_outfile_${CHIP}_femNPC.tab.gz --refPointLabel "peak center" --colorMap PuBuGn -T "${i} Clusters: ${CHIP}merge_femNPC signal" --outFileName ${outfolder}/linkagemethod_complete/heatmaps/${i}Clusters_merge_femNPC_${CHIP}.png --whatToShow "heatmap and colorbar" &	

done
done

#merge the ES and NPC's heatmaps.. inside folder linkagemethod_ward/heatmaps
#A) Clusterwise
for i in 7 8 9 10
do
montage ${i}Clusters_merge_femES_MCRS1.png ${i}Clusters_merge_femES_MOF.png ${i}Clusters_merge_femES_MSL2.png ${i}Clusters_merge_femES_NSL3.png ${i}Clusters_merge_femES_MSL1.png ${i}Clusters_merge_femNPC_MOF.png ${i}Clusters_merge_femNPC_MSL2.png ${i}Clusters_merge_femNPC_NSL3.png  \
-tile x1 -geometry +0.1+0.1 ChIPclusters${i}_HMcollage_ES_NPC_new.png
done

#B) Final 
montage ChIPclusters7_HMcollage_ES_NPC.png ChIPclusters8_HMcollage_ES_NPC.png ChIPclusters9_HMcollage_ES_NPC.png ChIPclusters10_HMcollage_ES_NPC.png -tile 1x -geometry +0.1+0.1 ChIPclusters_linkageWARD_Cluster7-10_ALL_HMcollage_ES_NPC.png

