#!/bin/bash

set -e

BAM_FILE=$1
QUALITY=35
if [[ -n "${2}" ]] ; then
    QUALITY=$2
fi

MIN_PAIRS=2
if [[ -n "${3}" ]] ; then
    MIN_PAIRS=$3
fi

WORKING_DIR=`dirname $BAM_FILE`
FILE_NAME=`basename $BAM_FILE`
FILE_PREFIX=`echo $FILE_NAME | awk -F. '{ print $1 }'`

CONFIG_FILE=${FILE_PREFIX}_${QUALITY}_${MIN_PAIRS}.cfg

bam2cfg.pl $WORKING_DIR/$BAM_FILE > $CONFIG_FILE
BD_OUT_FILE=${FILE_PREFIX}_${QUALITY}_${MIN_PAIRS}.bd_out.txt
breakdancer_max -d ${FILE_PREFIX}.bdreads.bed -g ${FILE_PREFIX}_${QUALITY}_${MIN_PAIRS}.bd.bed $CONFIG_FILE -q $QUALITY -r $MIN_PAIRS > $BD_OUT_FILE

#cat $BD_OUT_FILE | awk '!/^#/ {OFS="\t"; print $1,$2,$5,$7":"$8,$9,"+"}' > ${FILE_PREFIX}.bd_regions.bed
cat $BD_OUT_FILE | awk '!/^#/ && $1 == $4 {OFS="\t"; print $1,$2,$5,$7":"$8,$9,"+"}' > ${FILE_PREFIX}_${QUALITY}_${MIN_PAIRS}.bd_regions.bed
cat $BD_OUT_FILE | awk '!/^#/ && $1 != $4 {OFS="\t"; print $1,$2,$2+1000,$7"("$1"-"$4"):"$8,$9,"+"}' >> ${FILE_PREFIX}_${QUALITY}_${MIN_PAIRS}.bd_regions.bed
cat $BD_OUT_FILE | awk '!/^#/ && $1 != $4 {OFS="\t"; print $4,$5,$5+1000,$7"("$1"-"$4"):"$8,$9,"+"}' >> ${FILE_PREFIX}_${QUALITY}_${MIN_PAIRS}.bd_regions.bed
