#!/bin/bash

BAM_PATH=$1
WORKING_DIR=`dirname $BAM_PATH`
BAM_FILENAME=`basename $BAM_PATH`
FILE_PREFIX=`echo $BAM_FILENAME | awk -F. '{ print $1 }'`

echo "java -Xmx4g -jar /g/whelanch/software/picard-tools-1.72/CleanSam.jar I=$BAM_PATH O=$WORKING_DIR/${FILE_PREFIX}_clean.bam VALIDATION_STRINGENCY=LENIENT"
java -Xmx4g -jar /g/whelanch/software/picard-tools-1.72/CleanSam.jar I=$BAM_PATH O=$WORKING_DIR/${FILE_PREFIX}_clean.bam VALIDATION_STRINGENCY=LENIENT