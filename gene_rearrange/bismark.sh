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

echo bismark --path_to_bowtie /g/whelanch/software/condor_run_wrappers -q \
echo        $FORMAT --directional -n $MISMATCHES -l $SEED_LEN -e $MAQERR --chunkmbs 512 --unmapped --ambiguous \
echo	$REF -1 $READS1 -2 $READS2 

bismark --path_to_bowtie /g/whelanch/software/condor_run_wrappers -q \
        $FORMAT --directional -n $MISMATCHES -l $SEED_LEN -e $MAQERR --chunkmbs 512 --unmapped --ambiguous \
	$REF -1 $READS1 -2 $READS2 
