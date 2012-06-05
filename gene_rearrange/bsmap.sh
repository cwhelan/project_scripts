#!/bin/bash

# $1 - first sequence file
# $2 - second sequence file
# $3 - ref
# $4 - output_directory

OUTFILE=`basename $1`.bsmap.bsp
UNPAIRED=`basename $1`.unpaired.bsp

bsmap -a $1 -b $2 -d $3  \
	-o $4/$OUTFILE -2 $4/$UNPAIRED -p 12 -w 10 


