#!/bin/bash

set -e
set -u

REFERENCE=$1
CFG_FILE=$2
OUTPUT_PREFIX=$3

pindel -f $REFERENCE -i $CFG_FILE -c ALL -o $OUTPUT_PREFIX