#!/bin/bash

REF=$1
INPUT=$2
OUTPUT=$3
IGNORE_SUFFIX=$4

IGNORE_FILE=$REF$IGNORE_SUFFIX

REF_FILE=`basename $1`
REF_DIR=`dirname $1`

# pash-3.0lx.exe -h ref.dnameth.fa -v myReads.fastq -G 1 -k 13 -n 21 -o output2.txt -s 30 -d 800 -S . -B
pash-3.0lx.exe -h $1 -v $2 -G 1 -k 12 -n 18 -o $3 -s 30 -d 800 -S /l2/users/whelanch/gene_rearrange/pash/tmp -B -L $IGNORE_FILE

