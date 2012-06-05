#!/bin/bash

INPUT_BED=$1
TARGET_REGIONS=$2
OUTPUT=$3

WORKING_DIR=`dirname $INPUT_BED`

BED_FILENAME=`basename $INPUT_BED`
BAM_PARAM=`echo $BED_FILENAME | awk -F. '$2 == "bam" { print "-abam" }'`

if [[ -z $BAM_PARAM ]] ; then
    BAM_PARAM="-a"
fi


echo "coverageBed -d $BAM_PARAM $INPUT_BED -b $TARGET_REGIONS | $GR_HOME/scripts/compute_average_coverage_from_coverageBedD.py > $WORKING_DIR/$OUTPUT"
cat $TARGET_REGIONS | cut -f 1-3 | coverageBed -d $BAM_PARAM $INPUT_BED -b stdin | $GR_HOME/scripts/compute_average_coverage_from_coverageBedD.py > $WORKING_DIR/$OUTPUT
#cat $TARGET_REGIONS | cut -f 1-3 | coverageBed -d $BAM_PARAM $INPUT_BED -b stdin > $WORKING_DIR/$OUTPUT
