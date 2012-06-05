#!/bin/bash

samtools view -uf 0x2 -F 0x0400 novoalign_tier1_rg_sorted_clean_rmz_mdup_nsort.bam | bamToBed -i stdin -bedpe | awk '{OFS="\t"; print $1,$2,$6,$7,$8}' | awk '$3 - $2 > 500' | sort -k1,1 -k2,2n > novoalign_tier1_rg_sorted_clean_rmz_mdup_spans_li.bed
