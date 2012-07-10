#!/bin/bash

set -e
set -u

PARENT_DIR=$1

find $PARENT_DIR -name novoalign_tier1_sort_clean_mdup.bam | \
    sed 's/^/INPUT=/' |
    xargs java -jar /g/whelanch/software/picard-tools-1.72/MergeSamFiles.jar OUTPUT=$PARENT_DIR/novoalign_tier1_sort_clean_mdup.bam USE_THREADING=true