#!/bin/bash
#samtools view -uF 2 /l2/users/whelanch/gene_rearrange/tier1_bwa.bam | \
#           bamToFastq -bam stdin \
#                      -fq1 /l2/users/whelanch/gene_rearrange/s.tier1.disc.1.fq \
#                      -fq2 /l2/users/whelanch/gene_rearrange/s.tier1.disc.2.fq
set -e
set -u


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
SAMPLE=$8
RG=$9

PLATFORM="ILLUMINA"


/g/whelanch/software/bin/novoalign -d $REF \
                         -c 1 -f $f1 $f2 \
                         -F ILMFQ -k -K calfile.txt -i MP $5,$6 \
			 -r $REPEAT_REPORT -oSAM $"@RG\tID:$RG\tPU:$PLATFORM\tLB:$LIBRARY_NAME\tSM:$SAMPLE" 2> error.txt | \
	/g/whelanch/software/bin/samtools view -Sb - > ./novoalign.bam 
if (! grep "# Done at" error.txt); 
then 
    cat error.txt > /dev/stderr
    exit 1
fi
cat error.txt > /dev/stderr
/g/whelanch/software/bin/samtools sort -n ./novoalign.bam novoalign_sorted
