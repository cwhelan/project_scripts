#!/bin/bash

echo "/g/whelanch/software/bin/bwa sampe $1 $2 $3 $4 $5 | /g/whelanch/software/bin/samtools view -Sb - > $6"
/g/whelanch/software/bin/bwa sampe $1 $2 $3 $4 $5 | /g/whelanch/software/bin/samtools view -Sb - > $6
