#!/bin/bash

set -e
set -u

BAM_FILE=$1
REF=$2

MAPQ=0
if [[ -n "${3}" ]] ; then
    MAPQ=$3
fi

MAD_CUTOFF=5
if [[ -n "${4}" ]] ; then
    MAD_CUTOFF=$4
fi


WORKING_DIR=`dirname $BAM_FILE`
FILE_NAME=`basename $BAM_FILE`
FILE_PREFIX=`echo $FILE_NAME | awk -F. '{ print $1 }'`

PE_OUTFILE=${FILE_PREFIX}.delly_q${MAPQ}_c${MAD_CUTOFF}_del.txt
SR_OUTFILE=${FILE_PREFIX}.delly_q${MAPQ}_c${MAD_CUTOFF}_br.txt

echo "runing DELLY"
date +%s
echo "delly -p -g $REF -q $MAPQ -s $MAD_CUTOFF $BAM_FILE -o $PE_OUTFILE -b $SR_OUTFILE"
delly -p -g $REF -q $MAPQ -s $MAD_CUTOFF $BAM_FILE -o $PE_OUTFILE -b $SR_OUTFILE
echo "done"
date +%s