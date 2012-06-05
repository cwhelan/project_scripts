#!/bin/bash

set -e
set -u

OUTPUT=$1
WORKING_DIR=$2
IN_FILE_NAMES=$3

#echo "in file names: $IN_FILE_NAMES"

# in case there are too many chunk files, process them in batches of 100

NUM_BATCHES=1
NUM_CHUNK_FILES=0
CURRENT_BATCH_SIZE=0
CURRENT_BATCH=""
for CHUNK_FILE in wd*/$IN_FILE_NAMES
do
    #echo "chunk file $NUM_CHUNK_FILES $CHUNK_FILE"
    NUM_CHUNK_FILES=`echo "$NUM_CHUNK_FILES + 1" | bc`
    CURRENT_BATCH="$CURRENT_BATCH $CHUNK_FILE"
    CURRENT_BATCH_SIZE=`echo "$CURRENT_BATCH_SIZE + 1" | bc`
    if [ $NUM_CHUNK_FILES -eq 100 ]
    then
	echo "samtools merge -n tmp.${OUTPUT}.$NUM_BATCHES $CURRENT_BATCH"
	samtools merge -n tmp.${OUTPUT}.$NUM_BATCHES $CURRENT_BATCH
	CURRENT_BATCH=""
	NUM_CHUNK_FILES=0
	NUM_BATCHES=`echo "$NUM_BATCHES + 1" | bc`	
    fi
done

echo "samtools merge -n tmp.${OUTPUT}.$NUM_BATCHES $CURRENT_BATCH"
samtools merge -n tmp.${OUTPUT}.$NUM_BATCHES $CURRENT_BATCH

if [[ $NUM_BATCHES -gt 1 ]]
then
    echo "samtools merge -n $OUTPUT.bam $WORKING_DIR/tmp.$OUTPUT.*"
    samtools merge -n $OUTPUT.bam $WORKING_DIR/tmp.$OUTPUT.*
else
    echo "cp tmp.${OUTPUT}.1 ${OUTPUT}.bam"
    cp tmp.${OUTPUT}.1 ${OUTPUT}.bam
fi

rm $WORKING_DIR/tmp.$OUTPUT.*

#echo "samtools merge -n $OUTPUT $WORKING_DIR/wd*/$IN_FILE_NAMES" 1>&2
#samtools merge -n ${OUTPUT}.bam $WORKING_DIR/wd*/$IN_FILE_NAMES
samtools sort -m 5000000000 ${OUTPUT}.bam ${OUTPUT}_sort

