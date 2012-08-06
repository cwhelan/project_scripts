#!/bin/bash

WORKING_DIR=$1

rm -rf $WORKING_DIR/wd*
rm novoalign_tier1.bam
rm novoalign_tier1_sort.bam
rm novoalign_tier1_sort_clean.bam

rm novoalign_tier1_sort_clean_rmdup_long_inserts.txt
rm novoalign_tier1_sort_clean_rmdup_short_inserts.txt

gzip *.fq
gzip *.txt