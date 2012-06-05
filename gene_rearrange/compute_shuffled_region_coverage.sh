#!/bin/bash

INPUT_BED=$1
REF=$2
TARGET_REGIONS=$3
EXCLUDE_REGIONS=$4
OUTPUT=$5

WORKING_DIR=`dirname $INPUT_BED`

sleep $(($RANDOM % 600))

shuffleBed -i $INPUT_BED -g $REF.fai -excl $EXCLUDE_REGIONS -chrom | coverageBed -d -a stdin -b $TARGET_REGIONS | $GR_HOME/scripts/compute_average_coverage_from_coverageBedD.py > $WORKING_DIR/$OUTPUT
