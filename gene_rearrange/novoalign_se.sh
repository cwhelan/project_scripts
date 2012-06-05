#!/bin/bash

set -e 
set -u

f1=$1
REF=$2
WORKING_DIR=$3
REPEAT_REPORT=$4
if [[ -z "$4" ]]; then
    REPEAT_REPORT="Random"
fi
LIBRARY_NAME=$5

RG="RG1"
PLATFORM="ILLUMINA"
LIBRARY=$LIBRARY_NAME

F_PREFIX=`basename $f1 | awk -F. '{ print $1 }'`

/g/whelanch/software/bin/novoalign -d $REF \
                         -f $f1  \
                         -F ILMFQ -r $REPEAT_REPORT -oSAM $"@RG\tID:$RG\tPU:$PLATFORM\tLB:$LIBRARY\tSM:$LIBRARY" | \
	/g/whelanch/software/bin/samtools view -Sb - > $WORKING_DIR/$F_PREFIX.bam 
/g/whelanch/software/bin/samtools sort -n $WORKING_DIR/$F_PREFIX.bam $WORKING_DIR/${F_PREFIX}_sorted
