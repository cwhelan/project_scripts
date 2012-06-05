#!/bin/bash
#samtools view -uF 2 /l2/users/whelanch/gene_rearrange/tier1_bwa.bam | \
#           bamToFastq -bam stdin \
#                      -fq1 /l2/users/whelanch/gene_rearrange/s.tier1.disc.1.fq \
#                      -fq2 /l2/users/whelanch/gene_rearrange/s.tier1.disc.2.fq
f1=$1
f2=$2
REF=$3

RG="RG1"
PLATFORM="ILLUMINA"
LIBRARY="LC_1"

REPEAT_REPORT="Random"

/g/whelanch/software/bin/novoalign -d $REF \
                         -c 1 -f $f1 $f2 \
                         -F ILMFQ -k -i MP 5000,800 75,35 -a -r $REPEAT_REPORT -oSAM $"@RG\tID:$RG\tPU:$PLATFORM\tLB:$LIBRARY\tSM:$LIBRARY" | \
	/g/whelanch/software/bin/samtools view -Sb - > ./novoalign.bam 
/g/whelanch/software/bin/samtools sort -n ./novoalign.bam novoalign_sorted
