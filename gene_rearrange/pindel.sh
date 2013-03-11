#!/bin/bash

set -e
set -u

REFERENCE=$1
CFG_FILE=$2
OUTPUT_PREFIX=$3
ADDITIONAL_PINDEL_PARAMS=$4

echo "running pindel"
echo "pindel $ADDITIONAL_PINDEL_PARAMS -f $REFERENCE -i $CFG_FILE -c ALL -o $OUTPUT_PREFIX"
date +%s
pindel $ADDITIONAL_PINDEL_PARAMS -f $REFERENCE -i $CFG_FILE -c ALL -o $OUTPUT_PREFIX
echo "done"
date +%s