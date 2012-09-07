#!/bin/bash

set -e
set -u

PARENT_DIR=$1

find $PARENT_DIR -name novoalign_tier1_sort_clean_mdup.bam | \
    sed 's/^/INPUT=/' |
    xargs echo "java -Xmx4096m -jar /g/whelanch/software/picard-tools-1.72/MergeSamFiles.jar  OUTPUT=$PARENT_DIR/novoalign_tier1_sort_clean_mdup.bam USE_THREADING=true"

find $PARENT_DIR -name novoalign_tier1_sort_clean_mdup.bam | \
    sed 's/^/INPUT=/' |
    xargs java -Xmx4096m -jar /g/whelanch/software/picard-tools-1.72/MergeSamFiles.jar  OUTPUT=$PARENT_DIR/novoalign_tier1_sort_clean_mdup.bam USE_THREADING=true

find $PARENT_DIR -name novoalign_tier1_sort_clean_rmdup.bam | \
    sed 's/^/INPUT=/' | 
    xargs echo "java -Xmx4096m -jar /g/whelanch/software/picard-tools-1.72/MergeSamFiles.jar  OUTPUT=$PARENT_DIR/novoalign_tier1_sort_clean_rmdup.bam USE_THREADING=true"

find $PARENT_DIR -name novoalign_tier1_sort_clean_rmdup.bam | \
    sed 's/^/INPUT=/' | 
    xargs java -Xmx4096m -jar /g/whelanch/software/picard-tools-1.72/MergeSamFiles.jar OUTPUT=$PARENT_DIR/novoalign_tier1_sort_clean_rmdup.bam USE_THREADING=true