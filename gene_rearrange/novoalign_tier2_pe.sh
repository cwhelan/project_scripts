#!/bin/bash

set -e
set -u

f1=$1
f2=$2
REF=$3

RG="RG1"
PLATFORM="ILLUMINA"
LIBRARY="LC_1"

TARGET_ISIZE=$4
ISIZE_SD=$5
THRESHOLD=$6

CALFILE_PARAM=""
if [ -e ../calfile.txt ]
then
    CALFILE_PARAM="-k ../calfile.txt "
fi

/g/whelanch/software/bin/novoalign -d $REF \
                         -c 1 -f $f1 $f2 \
			 -a GATCGGAAGAGCGGTTCAGCA GATCGGAAGAGCGTCGTGTAGGGA \
                         $CALFILE_PARAM -i PE $TARGET_ISIZE,$ISIZE_SD -a -r Ex 1100 -t $THRESHOLD \
                         -oSAM $"@RG\tID:$RG\tPU:$PLATFORM\tLB:$LIBRARY\tSM:$LIBRARY" | \
	/g/whelanch/software/bin/samtools view -Sb - > ./novoalign.bam 
/g/whelanch/software/bin/samtools sort -n ./novoalign.bam novoalign_sorted
