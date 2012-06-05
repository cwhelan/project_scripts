#!/bin/bash

READ_FILE1=$1
READ_FILE2=$2

SCRIPTS=/l2/users/whelanch/gene_rearrange/scripts

$SCRIPTS/cutadapt.sh $READ_FILE1 $READ_FILE1.cut.fq
$SCRIPTS/cutadapt.sh $READ_FILE2 $READ_FILE2.cut.fq

$SCRIPTS/discard_short_reads.py $READ_FILE1.cut.fq $READ_FILE2.cut.fq $READ_FILE1.trim.fq $READ_FILE2.trim.fq 48 101
rm $READ_FILE1.cut.fq
rm $READ_FILE2.cut.fq
