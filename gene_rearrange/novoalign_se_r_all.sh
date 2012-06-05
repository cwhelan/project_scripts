#!/bin/bash

f1=$1
REF=$2
WORKING_DIR=$3

RG="RG1"
PLATFORM="ILLUMINA"
LIBRARY="LC_1"

F_PREFIX=`basename $f1 | awk -F. '{ print $1 }'`

/g/whelanch/software/bin/novoalign -d $REF \
                         -f $f1  \
                         -F ILMFQ -r All -oSAM $"@RG\tID:$RG\tPU:$PLATFORM\tLB:$LIBRARY\tSM:$LIBRARY" > $WORKING_DIR/${F_PREFIX}_r_all.sam
#	/g/whelanch/software/bin/samtools view -Sb - > $WORKING_DIR/$F_PREFIX_r_all.bam 
