#!/bin/bash

set -e
set -u

REF=$1

for i in wd* ; 
do 
samtools sort $i/novoalign.CT.bam $i/novoalign.CT.sort;
samtools sort $i/novoalign.GA.bam $i/novoalign.GA.sort; 
done

$GR_HOME/scripts/merge_bams.sh novoalign.CT . novoalign.CT.sort.bam
$GR_HOME/scripts/merge_bams.sh novoalign.GA . novoalign.GA.sort.bam

echo "samtools mpileup -BC 0 -q 30 -f $REF novoalign.CT_sort.bam novoalign.GA_sort.bam > novoalign.mpileup.txt"
samtools mpileup -BC 0 -q 30 -f $REF novoalign.CT_sort.bam novoalign.GA_sort.bam > novoalign.mpileup.txt

# homogenous
novomethyl -h <novoalign.mpileup.txt > novoalign.Cmethyl.hom.bed

# heterogenous
novomethyl -h -% <novoalign.mpileup.txt > novoalign.Cmethyl.het.bed

# consensus
novomethyl -h -o Consensus <novoalign.mpileup.txt > novoalign.Cmethyl.cons.txt
