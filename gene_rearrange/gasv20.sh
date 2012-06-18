#!/bin/bash

set -e
set -u

BAM_FILE=$1

GASV_HOME=/g/whelanch/software/gasv
GASVPRO_HOME=/g/whelanch/software/GASVProSrc

FILE_PREFIX=`echo $BAM_FILE | awk -F. '{ print $1 }'`

echo "java -Xms512m -Xmx2048m -jar ${GASV_HOME}/bin/BAMToGASV.jar $BAM_FILE -LIBRARY_SEPARATED all"
java -Xms512m -Xmx2048m -jar ${GASV_HOME}/bin/BAMToGASV.jar $BAM_FILE -LIBRARY_SEPARATED all

echo "java -Xms512m -Xmx2048m -jar ${GASV_HOME}/bin/GASV.jar --minClusterSize 3 --output 2"
java -Xms512m -Xmx2048m -jar ${GASV_HOME}/bin/GASV.jar --minClusterSize 3 --output 2

#echo "Pruning Redundant/Overlapping Predictions...."
#${GASVPRO_HOME}/scripts/pruneClusters.pl ${FILE_PREFIX}_${LIBRARY_NAME}.deletion.clusters ESP > prune_clusters.out



