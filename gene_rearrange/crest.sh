#!/bin/bash

set -e
set -u

BAM_FILE=$1
REF=$2
CHROM=$3

/g/whelanch/software/bin/extractSClip.pl -i $BAM_FILE --ref $REF -r $CHROM
