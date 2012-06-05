#!/bin/bash

set -e
set -u

BAM_FILE=$1

WORKING_DIR=`dirname $BAM_FILE`
FILE_NAME=`basename $BAM_FILE`
FILE_PREFIX=`echo $FILE_NAME | awk -F. '{ print $1 }'`

LI_BAM_FILE=${FILE_PREFIX}_li.bam
LI_CONFIG_FILE=${FILE_PREFIX}_li.cfg

# separate long inserts
samtools view -h $BAM_FILE | grep -v 'ZO:Z:+-' | samtools view -Sb - > $WORKING_DIR/$LI_BAM_FILE

bam2cfg.pl $WORKING_DIR/$LI_BAM_FILE > $LI_CONFIG_FILE
BD_OUT_FILE=${FILE_PREFIX}_li.bd_out.txt
breakdancer_max -a -d ${FILE_PREFIX}_li.bdreads.bed -g ${FILE_PREFIX}_li.bd.bed -l $LI_CONFIG_FILE > $BD_OUT_FILE

cat $BD_OUT_FILE | awk '!/^#/ && $1 == $4 {OFS="\t"; print $1,$2,$5,$7":"$8,$9,"+"}' > ${FILE_PREFIX}_li.bd_regions.bed
cat $BD_OUT_FILE | awk '!/^#/ && $1 != $4 {OFS="\t"; print $1,$2,$2+1000,$7"("$1"-"$4"):"$8,$9,"+"}' >> ${FILE_PREFIX}.bd_regions.bed
cat $BD_OUT_FILE | awk '!/^#/ && $1 != $4 {OFS="\t"; print $4,$5,$5+1000,$7"("$1"-"$4"):"$8,$9,"+"}' >> ${FILE_PREFIX}.bd_regions.bed
