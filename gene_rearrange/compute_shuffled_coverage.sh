#!/bin/bash

INPUT_BED=$1
REF=$2
TARGET_REGIONS=$3
OUTPUT=$4

shuffleBed -i novoalign_tier1_rg_sorted_clean_rmz_mdup_spans.bed -g /l2/users/whelanch/genome_refs/10KG/hg19/hg19.fa.fai -excl /l2/users/whelanch/genome_refs/10KG/hg19/gaps.bed -chrom | coverageBed -a stdin -b /l2/users/whelanch/genome_refs/10KG/hg19/repeatmasker_merged.bed > repeat_coverage_shuffled.bed