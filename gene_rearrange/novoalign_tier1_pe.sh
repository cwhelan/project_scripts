#!/bin/bash
#samtools view -uF 2 /l2/users/whelanch/gene_rearrange/tier1_bwa.bam | \
#           bamToFastq -bam stdin \
#                      -fq1 /l2/users/whelanch/gene_rearrange/s.tier1.disc.1.fq \
#                      -fq2 /l2/users/whelanch/gene_rearrange/s.tier1.disc.2.fq
set -e
set -u

source 'novo_params'

echo "sourced novo_params"
echo "FORMAT = $FORMAT"

f1=$1
f2=$2
REF=$3

REPEAT_REPORT=$4
if [[ -z "$4" ]]; then
    REPEAT_REPORT="Random"
fi

TARGET_ISIZE=$5
ISIZE_SD=$6
LIBRARY_NAME=$7

RG="RG1"
PLATFORM="ILLUMINA"
LIBRARY=${LIBRARY_NAME}


/g/whelanch/software/bin/novoalign -d $REF \
                         -c 1 -f $f1 $f2 \
                         -F $FORMAT -k -K calfile.txt -i PE $5,$6 -a -r $REPEAT_REPORT -oSAM $"@RG\tID:$RG\tPU:$PLATFORM\tLB:$LIBRARY\tSM:$LIBRARY" | \
	/g/whelanch/software/bin/samtools view -Sb - > ./novoalign.bam 
/g/whelanch/software/bin/samtools sort -n ./novoalign.bam novoalign_sorted
