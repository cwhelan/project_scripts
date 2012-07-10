#!/bin/bash

WORKING_DIR=$1

rm -rf $WORKING_DIR/wd*
rm novoalign_tier1.bam
rm novoalign_tier1_sort.bam
rm novoalign_tier1_sort_clean.bam

gzip *.fq
gzip *.txt