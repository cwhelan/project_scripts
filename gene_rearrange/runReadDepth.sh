#!/bin/sh

set -e
set -u

BAMFILE=$1
READ_LENGTH=$2

[ -f params ]

mkdir annotations
mkdir annotations/gcWinds
ln -s $L2H/readDepth/annotations/hg19/gcWinds/readLength$READ_LENGTH/* annotations/gcWinds/

mkdir annotations/mapability
ln -s $L2H/readDepth/annotations/hg19/mapability/readLength$READ_LENGTH/* annotations/mapability/

ln -s $L2H/readDepth/annotations/hg19/entrypoints annotations/

mkdir output

mkdir reads

bamToBed -i $BAMFILE | translateNCBIChromNamesToUCSC -bed | awk '{close(f);f=$1}{print > "reads/"f".bed"}'

Rscript $GR_HOME/scripts/readDepth.R

echo "track name=COPY_NUMBER_SEGMENTS" > seg_track_line
cat output/segs.dat | awk '{OFS="\t"; print $1,$2,$3,$5}' | cat seg_track_line - > output/copy_number_segments.bed

echo "track name=ALTERED_COPY_NUMBER" > alt_track_line
cat output/alts.dat | awk '{OFS="\t"; print $1,$2,$3,$5}' | cat alt_track_line - > output/altered_copy_number.bed