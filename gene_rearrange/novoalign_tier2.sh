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

/g/whelanch/software/bin/novoalign -d $REF \
                         -c 1 -f $f1 $f2 \
                         -k ../calfile.txt -i MP $TARGET_ISIZE,$ISIZE_SD 150,50 \
			 -a GATCGGAAGAGCGGTTCAGCA GATCGGAAGAGCGTCGTGTAGGGA \
			 -r Ex 1100 -t $THRESHOLD -oSAM $"@RG\tID:$RG\tPU:$PLATFORM\tLB:$LIBRARY\tSM:$LIBRARY" 2> error.txt | \
	/g/whelanch/software/bin/samtools view -Sb - > ./novoalign.bam
if (! grep "# Done at" error.txt); 
then 
    cat error.txt > /dev/stderr
    exit 1
fi
cat error.txt > /dev/stderr
/g/whelanch/software/bin/samtools sort -n ./novoalign.bam novoalign_sorted
