#!/bin/bash

#expects a name-sorted bed
# found here: http://groups.google.com/group/bedtools-discuss/web/community-usage-examples
INPUT_BAM=$1
REF=$2

samtools view -uf 0x2 $INPUT_BAM | bamToBed -i stdin -bedpe | awk '{OFS="\t"; print $1,$2,$6,$7,$8}' | sort -k1,1 -k2,2n > $INPUT_BAM.bed
genomeCoverageBed -i stdin -g $REF > $INPUT_BAM.spanning_coverage.txt
