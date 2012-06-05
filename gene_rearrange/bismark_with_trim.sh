#!/bin/bash

set -e
set -u

REF=$1
READS1=$2
READS2=$3
MISMATCHES=$4
SEED_LEN=$5
MAQERR=$6
FORMAT=$7

PREFIX=`echo $READS1 | awk -F.gz '{ print $1 }'`

/l2/users/whelanch/gene_rearrange/scripts/trim_adapters_pe.sh $READS1 $READS2

echo bismark --path_to_bowtie /g/whelanch/software/condor_run_wrappers -q \
echo        $FORMAT --directional -n $MISMATCHES -l $SEED_LEN -e $MAQERR --chunkmbs 512 --unmapped --ambiguous \
echo	$REF -1 $READS1.trim.fq -2 $READS2.trim.fq

bismark --path_to_bowtie /g/whelanch/software/condor_run_wrappers -q \
        $FORMAT --directional -n $MISMATCHES -l $SEED_LEN -e $MAQERR --chunkmbs 512 --unmapped --ambiguous \
	$REF -1 $READS1.trim.fq -2 $READS2.trim.fq
