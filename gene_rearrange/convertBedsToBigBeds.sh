#!/bin/bash

set -e
set -u

for f in *.bed
do
    s=`du $f | cut -f1`
    if [ "$s" -ge 50000 ];
    then
	prefix=`echo $f | awk -F. '{ print $1 }'`
	echo bedToBigBed $prefix.bed /l2/users/whelanch/genome_refs/10KG/hg19/hg19.fa.fai $prefix.bb
	bedToBigBed $prefix.bed /l2/users/whelanch/genome_refs/10KG/hg19/hg19.fa.fai $prefix.bb
	rm $prefix.bed
    else
	echo $f is too small to big-bedify.
    fi
done