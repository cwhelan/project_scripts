#!/bin/bash

BAM_PATH=$1
WORKING_DIR=`dirname $BAM_PATH`
BAM_FILENAME=`basename $BAM_PATH`
FILE_PREFIX=`echo $BAM_FILENAME | awk -F. '{ print $1 }'`

echo "java -Xmx4g -jar /g/whelanch/software/picard-tools-1.40/MarkDuplicates.jar I=$BAM_PATH O=$WORKING_DIR/${FILE_PREFIX}_mdup.bam METRICS_FILE=$WORKING_DIR/${FILE_PREFIX}_dup_metrics.txt ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT"
java -Xmx8g -jar /g/whelanch/software/picard-tools-1.40/MarkDuplicates.jar I=$BAM_PATH O=$WORKING_DIR/${FILE_PREFIX}_mdup.bam METRICS_FILE=$WORKING_DIR/${FILE_PREFIX}_dup_metrics.txt ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT

echo "java -Xmx4g -jar /g/whelanch/software/picard-tools-1.40/MarkDuplicates.jar I=$BAM_PATH O=$WORKING_DIR/${FILE_PREFIX}_mdup.bam METRICS_FILE=$WORKING_DIR/${FILE_PREFIX}_dup_metrics.txt ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true"
java -Xmx8g -jar /g/whelanch/software/picard-tools-1.40/MarkDuplicates.jar I=$BAM_PATH O=$WORKING_DIR/${FILE_PREFIX}_rmdup.bam METRICS_FILE=$WORKING_DIR/${FILE_PREFIX}_dup_metrics.txt ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true

samtools sort -n $WORKING_DIR/${FILE_PREFIX}_mdup.bam $WORKING_DIR/${FILE_PREFIX}_mdup_nsort