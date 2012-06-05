#!/bin/bash

set -e
set -u

INPUT_BAM=$1

WORKING_DIR=`dirname $INPUT_BAM`

BAM_FILENAME=`basename $INPUT_BAM`

echo "samtools view -uF 0x2 $INPUT_BAM | \
           bamToFastq -bam stdin \
                      -fq1 $WORKING_DIR/${BAM_FILENAME}_disc_1.fq \
                      -fq2 $WORKING_DIR/${BAM_FILENAME}_disc_2.fq"

samtools view -uF 0x2 $INPUT_BAM | \
           bamToFastq -bam stdin \
                      -fq1 $WORKING_DIR/${BAM_FILENAME}_disc_1.fq \
                      -fq2 $WORKING_DIR/${BAM_FILENAME}_disc_2.fq

