#!/bin/bash

set -e
set -u

BAM_FILE=$1
FILE_PREFIX=`echo $BAM_FILE | awk -F. '{ print $1 }'`

samtools flagstat $BAM_FILE > ${FILE_PREFIX}.flagstat.txt