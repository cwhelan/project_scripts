#!/bin/bash

set -e
set -u

BAMFILE=$1
FILE_PREFIX=`echo $BAMFILE | awk -F. '{ print $1 }'`

samtools view -F 0x40C $BAMFILE | awk '/ZO\:Z\:\-\+/ {print $9}'  > ${FILE_PREFIX}_long_inserts.txt
samtools view -F 0x40C $BAMFILE | awk '/ZO\:Z\:\+\-/ {print $9}' > ${FILE_PREFIX}_short_inserts.txt

echo "$GR_HOME/scripts/calculateInsertSizes.R ${FILE_PREFIX}_short_inserts.txt ${FILE_PREFIX}_long_inserts.txt > ${FILE_PREFIX}_insert_stats.txt"
$GR_HOME/scripts/calculateInsertSizes.R ${FILE_PREFIX}_short_inserts.txt ${FILE_PREFIX}_long_inserts.txt > ${FILE_PREFIX}_insert_stats.txt