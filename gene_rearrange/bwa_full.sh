#!/bin/bash

set -e
set -u

READ_FILE1=$1
READ_FILE2=$2
REFERENCE=$3
OUTFILE_PREFIX=$4
RGID=$5
LIBRARY=$6
SAMPLE=$7
ADDITIONAL_BWA_ALN_PARAMS=$8

READ1_BASE=`basename $READ_FILE1`
READ2_BASE=`basename $READ_FILE2`

echo "/g/whelanch/software/bin/bwa aln -t 8 $ADDITIONAL_BWA_ALN_PARAMS $REFERENCE $READ_FILE1 > $READ1_BASE.sai"
/g/whelanch/software/bin/bwa aln -t 8 $ADDITIONAL_BWA_ALN_PARAMS $REFERENCE $READ_FILE1 > $READ1_BASE.sai
echo "/g/whelanch/software/bin/bwa aln -t 8 $ADDITIONAL_BWA_ALN_PARAMS $REFERENCE $READ_FILE2 > $READ2_BASE.sai"
/g/whelanch/software/bin/bwa aln -t 8 $ADDITIONAL_BWA_ALN_PARAMS $REFERENCE $READ_FILE2 > $READ2_BASE.sai


echo "/g/whelanch/software/bin/bwa sampe -r ""'"@RG\tID:${RGID}\tLB:${LIBRARY}\tSM:${SAMPLE}"'"" $REFERENCE $READ1_BASE.sai $READ2_BASE.sai $READ_FILE1 $READ_FILE2 | samtools view -Sb - > $OUTFILE_PREFIX.bam"
/g/whelanch/software/bin/bwa sampe -r ""'"@RG\tID:${RGID}\tLB:${LIBRARY}\tSM:${SAMPLE}"'"" $REFERENCE $READ1_BASE.sai $READ2_BASE.sai $READ_FILE1 $READ_FILE2 | samtools view -Sb - > $OUTFILE_PREFIX.bam

echo samtools sort -m10000000000 $OUTFILE_PREFIX.bam ${OUTFILE_PREFIX}_sort
samtools sort -m10000000000 $OUTFILE_PREFIX.bam ${OUTFILE_PREFIX}_sort

echo samtools index ${OUTFILE_PREFIX}_sort.bam
samtools index ${OUTFILE_PREFIX}_sort.bam
