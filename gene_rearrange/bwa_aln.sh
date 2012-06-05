#!/bin/bash

echo "/g/whelanch/software/bin/bwa aln -t $5 -n $1 $2 $3 > $4"
/g/whelanch/software/bin/bwa aln -t $5 -n $1 $2 $3 > $4
