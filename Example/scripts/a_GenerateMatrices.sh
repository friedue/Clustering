#!/bin/sh
set -e # to make the script fail if a variable is not set
#######################################
# This script uses computeMatrix of the deepTools suite to calculate a 
# matrix of values (e.g. ChIP signals) for a  given set of genomic 
# regions. The matrix should be calculated for each factor. 
#######################################
## ${1} Path to signal files (bigWigs)
## ${2} Path to region files (bed)
## ${3} Folder for output
## ${4} Suffix for output files, e.g. X_ESC
######################################
##!!!!!!!!!!!!!!!
## make additional parameters of computeMatrix variable,
## e.g. the bin size (-bs), the region length (-m) etc. 

/package/deepTools-1.5.9.1/bin/computeMatrix scale-regions\
     -S ${1} \
     -R ${2} \
     -bs 50 -m 1500 --sortRegions no --missingDataAsZero \
     --outFileNameMatrix ${3}/Matrix_${4}.tab \
      --outFileName Matrix_${4}.gz
     # --outFileNameMatrix will always produce a gzipped file that is 
     # needed for heatmapper - actually, we don't need this output at
     # this point, but it is a non-optional parameter of computeMatrix
     # the output needed for the next steps is the one obtained with
     # --outFileNameMatrix
