#!/bin/bash

set -e
set -u

BAM_FILE=$1

GASV_HOME=/g/whelanch/software/gasv
GASVPRO_HOME=/g/whelanch/software/GASVProSrc

FILE_PREFIX=`echo $BAM_FILE | awk -F. '{ print $1 }'`

echo "java -Xms512m -Xmx4096m -jar ${GASV_HOME}/bin/BAMToGASV.jar $BAM_FILE"
java -Xms512m -Xmx4096m -jar ${GASV_HOME}/bin/BAMToGASV.jar $BAM_FILE

echo "java -Xms512m -Xmx4096m -jar ${GASV_HOME}/bin/GASV.jar --minClusterSize 3 --output reads ${BAM_FILE}.gasv.in"
java -Xms512m -Xmx4096m -jar ${GASV_HOME}/bin/GASV.jar --minClusterSize 3 --output reads --batch ${BAM_FILE}.gasv.in

#echo "Pruning Redundant/Overlapping Predictions...."
${GASVPRO_HOME}/scripts/pruneClusters.pl ${BAM_FILE}.gasv.in.clusters ESP > ${BAM_FILE}.gasv.in.clusters.pruned.clusters



