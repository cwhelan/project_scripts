#!/bin/bash

set -e
set -u

BAM_FILE=$1
FILE_PREFIX=`echo $BAM_FILE | awk -F. '{ print $1 }'`

echo "samtools flagstat $BAM_FILE > ${FILE_PREFIX}.flagstat.txt"
samtools flagstat $BAM_FILE > ${FILE_PREFIX}.flagstat.txt