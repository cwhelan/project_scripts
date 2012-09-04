#!/bin/bash

set -e
set -u

READ_FILE1=$1
READ_FILE2=$2
REFERENCE=$3
OUTFILE_NAME=$4

READ1_BASE=`basename $READ_FILE1`
READ2_BASE=`basename $READ_FILE2`

echo /g/whelanch/software/bin/bwa aln -t 8 $REFERENCE $READ_FILE1 > $READ1_BASE.sai
/g/whelanch/software/bin/bwa aln -t 8 $REFERENCE $READ_FILE1 > $READ1_BASE.sai
echo /g/whelanch/software/bin/bwa aln -t 8 $REFERENCE $READ_FILE2 > $READ2_BASE.sai
/g/whelanch/software/bin/bwa aln -t 8 $REFERENCE $READ_FILE2 > $READ2_BASE.sai
echo /g/whelanch/software/bin/bwa sampe $REFERENCE $READ1_BASE.sai $READ2_BASE.sai $READ_FILE1 $READ_FILE2 | samtools view -Sb - > $OUTFILE_NAME
/g/whelanch/software/bin/bwa sampe $REFERENCE $READ1_BASE.sai $READ2_BASE.sai $READ_FILE1 $READ_FILE2 | samtools view -Sb - > $OUTFILE_NAME
