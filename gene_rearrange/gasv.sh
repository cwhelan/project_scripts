#!/bin/bash

set -e
set -u

BAM_FILE=$1
LIBRARY_NAME=$2

GASV_HOME=/g/whelanch/software/gasv
GASVPRO_HOME=/g/whelanch/software/GASVProSrc

FILE_PREFIX=`echo $BAM_FILE | awk -F. '{ print $1 }'`

perl $GASV_HOME/bin/BAM_preprocessor.pl $BAM_FILE

INFO_FILE=$FILE_PREFIX.info 
# echo "cat $INFO_FILE | grep ${LIBRARY_NAME} | awk '{print $ 2}'"
#echo `cat $INFO_FILE | grep ${LIBRARY_NAME}`
#echo `cat $INFO_FILE | grep ${LIBRARY_NAME} | awk '{print $ 2}'`
L_MIN=`cat $INFO_FILE | grep -e ${LIBRARY_NAME} | awk '{print $ 2}'`
L_MAX=`cat $INFO_FILE | grep -e ${LIBRARY_NAME} | awk '{print $ 3}'`

echo "L_MIN: $L_MIN"
echo "L_MAX: $L_MAX"

java -jar ${GASV_HOME}/lib/gasv.jar --cluster --lmin $L_MIN --lmax $L_MAX  ${FILE_PREFIX}_${LIBRARY_NAME}.deletion

echo "Pruning Redundant/Overlapping Predictions...."
${GASVPRO_HOME}/scripts/pruneClusters.pl ${FILE_PREFIX}_${LIBRARY_NAME}.deletion.clusters ESP > prune_clusters.out



