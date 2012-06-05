#!/bin/bash

set -e
set -u 

SORT_BAM_PATH=$1
NSORT_BAM_PATH=$2
REF_INDEX_PATH=$3

WORKING_DIR=`dirname $SORT_BAM_PATH`
SORT_BAM_FILENAME=`basename $SORT_BAM_PATH`
SORT_FILE_PREFIX=`echo $SORT_BAM_FILENAME | awk -F. '{ print $1 }'`

NSORT_BAM_FILENAME=`basename $NSORT_BAM_PATH`
NSORT_FILE_PREFIX=`echo $NSORT_BAM_FILENAME | awk -F. '{ print $1 }'`

# echo "Coord. sorted bam file path: $SORT_BAM_PATH"
# echo "Name sorted bam file path: $NSORT_BAM_PATH"
# echo "Coord. sort file prefix: $SORT_FILE_PREFIX"
# echo "Name sort file prefix: $NSORT_FILE_PREFIX"

# # Find reads that are part of proper pairs (-f 0x2) and are not PCR duplicates (-F 0x400)

echo "samtools view -u -f 0x2 -F 0x400 $SORT_BAM_PATH | genomeCoverageBed -ibam stdin -g $REF_INDEX_PATH > $WORKING_DIR/${SORT_FILE_PREFIX}_read_coverage_hist.txt"
samtools view -u -f 0x2 -F 0x400 $SORT_BAM_PATH | genomeCoverageBed -ibam stdin -g $REF_INDEX_PATH > $WORKING_DIR/${SORT_FILE_PREFIX}_read_coverage_hist.txt

echo "samtools view -u -f 0x2 -F 0x400 $SORT_BAM_PATH | genomeCoverageBed -ibam stdin -g $REF_INDEX_PATH -bg > $WORKING_DIR/${SORT_FILE_PREFIX}_read_coverage.bed"
samtools view -u -f 0x2 -F 0x400 $SORT_BAM_PATH | genomeCoverageBed -ibam stdin -g $REF_INDEX_PATH -bg > $WORKING_DIR/${SORT_FILE_PREFIX}_read_coverage.bed

echo "samtools view -u -f 0x2 -f 0x400 $SORT_BAM_PATH | genomeCoverageBed -ibam stdin -g $REF_INDEX_PATH -bg > $WORKING_DIR/${SORT_FILE_PREFIX}_dup_coverage.bed"
samtools view -u -f 0x2 -f 0x400 $SORT_BAM_PATH | genomeCoverageBed -ibam stdin -g $REF_INDEX_PATH -bg > $WORKING_DIR/${SORT_FILE_PREFIX}_dup_coverage.bed

# Spanning coverage taken from http://groups.google.com/group/bedtools-discuss/web/community-usage-examples

# 1. Sort the BAM files by read ID so that step 2 can correctly compute the "span" of the
# paired-end fragment.
#echo "samtools sort -n -m 50000000000 $BAM_PATH $WORKING_DIR/${FILE_PREFIX}_nsort"
#samtools sort -n -m 50000000000 $BAM_PATH $WORKING_DIR/${FILE_PREFIX}_nsort

# 2. Create BED from the spanning coordinates and sort them by chrom for genomeCoverageBed.
# NOTE: Using -bedpe format from bamToBed so that I have the starts/ends for each mate.
# NOTE: Only computing coverage based upon "properly-paired" reads. (i.e., -f 0x2)
# NOTE: removing PCR duplicates from alignments (-F 0x400)

echo "samtools view -uf 0x2 -F 0x400 $NSORT_BAM_PATH | bamToBed -i stdin -bedpe | awk '{OFS="\t"; print \$1,\$2,\$6,\$7,\$8}' | sort -k1,1 -k2,2n > $WORKING_DIR/${NSORT_FILE_PREFIX}_spans.bed"
samtools view -uf 0x2 -F 0x400 $NSORT_BAM_PATH | bamToBed -i stdin -bedpe | awk '{OFS="\t"; print $1,$2,$6,$7,$8}' | sort -k1,1 -k2,2n > $WORKING_DIR/${NSORT_FILE_PREFIX}_spans.bed

# 3. Compute the ___spanning___ coverage histogram for the sample
echo "genomeCoverageBed -i $WORKING_DIR/${NSORT_FILE_PREFIX}_spans.bed -g $REF_INDEX_PATH > $WORKING_DIR/${NSORT_FILE_PREFIX}_span_coverage_hist.txt"
genomeCoverageBed -i $WORKING_DIR/${NSORT_FILE_PREFIX}_spans.bed -g $REF_INDEX_PATH > $WORKING_DIR/${NSORT_FILE_PREFIX}_span_coverage_hist.txt

echo "genomeCoverageBed -i $WORKING_DIR/${NSORT_FILE_PREFIX}_spans.bed -g $REF_INDEX_PATH -bg > $WORKING_DIR/${NSORT_FILE_PREFIX}_span_coverage.bed"
genomeCoverageBed -i $WORKING_DIR/${NSORT_FILE_PREFIX}_spans.bed -g $REF_INDEX_PATH -bg > $WORKING_DIR/${NSORT_FILE_PREFIX}_span_coverage.bed


# the same, but for uniquely aligned read coverage only

echo "samtools view -f 0x2 -F 0x400 $NSORT_BAM_PATH | grep -v 'ZS:Z:R' | samtools view -Sbu - | bamToBed -i stdin -bedpe | awk '{OFS="\t"; print \$1,\$2,\$6,\$7,\$8}' | sort -k1,1 -k2,2n > $WORKING_DIR/${NSORT_FILE_PREFIX}_uniq_spans.bed"
samtools view -h -f 0x2 -F 0x400 $NSORT_BAM_PATH | grep -v 'ZS:Z:R' | samtools view -Sbu - | bamToBed -i stdin -bedpe | awk '{OFS="\t"; print $1,$2,$6,$7,$8}' | sort -k1,1 -k2,2n > $WORKING_DIR/${NSORT_FILE_PREFIX}_uniq_spans.bed

# 3. Compute the ___spanning___ coverage histogram for the sample
echo "genomeCoverageBed -i $WORKING_DIR/${NSORT_FILE_PREFIX}_uniq_spans.bed -g $REF_INDEX_PATH > $WORKING_DIR/${NSORT_FILE_PREFIX}_uniq_span_coverage_hist.txt"
genomeCoverageBed -i $WORKING_DIR/${NSORT_FILE_PREFIX}_uniq_spans.bed -g $REF_INDEX_PATH > $WORKING_DIR/${NSORT_FILE_PREFIX}_uniq_span_coverage_hist.txt

echo "genomeCoverageBed -i $WORKING_DIR/${NSORT_FILE_PREFIX}_spans.bed -g $REF_INDEX_PATH -bg > $WORKING_DIR/${NSORT_FILE_PREFIX}_span_coverage.bed"
genomeCoverageBed -i $WORKING_DIR/${NSORT_FILE_PREFIX}_uniq_spans.bed -g $REF_INDEX_PATH -bg > $WORKING_DIR/${NSORT_FILE_PREFIX}_uniq_span_coverage.bed

