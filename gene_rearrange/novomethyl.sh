#!/bin/bash
set -e
set -u


f1=$1
f2=$2
REF=$3

TARGET_ISIZE=$4
ISIZE_SD=$5
LIBRARY_NAME=$6

RG="RG1"
PLATFORM="ILLUMINA"
LIBRARY=${LIBRARY_NAME}


/g/whelanch/software/bin/novoalign -b2 -d $REF.nbx -t130 -h120  \
                         -c 1 -f $f1 $f2 -a \
                         -F STDFQ -k -K calfile.txt -i PE $TARGET_ISIZE,$ISIZE_SD \
			 -oSAM $"@RG\tID:$RG\tPU:$PLATFORM\tLB:$LIBRARY\tSM:$LIBRARY" \
                         > novoalign.sam 2> novoalign.log
grep -e "^@\|ZB:Z:CT" novoalign.sam | samtools view -ubS - | samtools sort - novoalign.CT
grep -e "^@\|ZB:Z:GA" novoalign.sam | samtools view -ubS - | samtools sort - novoalign.GA
samtools mpileup -BC 0 -q 30 -f $REF novoalign.CT.bam novoalign.GA.bam > novoalign.mpileup.txt

# homogenous
novomethyl -h <novoalign.mpileup.txt > novoalign.Cmethyl.hom.bed

# heterogenous
novomethyl -h -% <novoalign.mpileup.txt > novoalign.Cmethyl.het.bed

# consensus
novomethyl -h -o Consensus <novoalign.mpileup.txt > novoalign.Cmethyl.cons.txt