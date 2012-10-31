#!/bin/bash

set -e
set -u

READGROUP=$1
LIBRARY=$2
SAMPLE=$3

cd $READGROUP
/l2/users/whelanch/project_scripts/gene_rearrange/build_novoalign_tier1.py `pwd` /l2/users/whelanch/genome_refs/10KG/hg19/hg19.fa.nix /l2/users/whelanch/genome_refs/10KG/hg19/hg19.fa.fai /l2/users/whelanch/1000genomes/pilot_data/data/$SAMPLE/sequence_read/${READGROUP}_1.recal.fastq.gz /l2/users/whelanch/1000genomes/pilot_data/data/$SAMPLE/sequence_read/${READGROUP}_2.recal.fastq.gz 200 30 Random $LIBRARY $SAMPLE $READGROUP PE STDFQ
mkdir /tmp/$READGROUP
cp novoalign_first_tier.dag /tmp/$READGROUP/
condor_submit_dag /tmp/$READGROUP/novoalign_first_tier.dag
cd ..
