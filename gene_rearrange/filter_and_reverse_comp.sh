#!/bin/bash
samtools view -h $1/$2 | /l2/users/whelanch/gene_rearrange/scripts/filter_short_inserts.py | samtools view -Sb - | \
           bamToFastq -bam stdin \
                      -fq1 $1/s_3_1.filter.txt \
                      -fq2 $1/s_3_2.filter.txt
/l2/users/whelanch/gene_rearrange/scripts/rev_comp.py $1/s_3_1.filter.txt $1/s_3_1.filter.rc.txt
/l2/users/whelanch/gene_rearrange/scripts/rev_comp.py $1/s_3_2.filter.txt $1/s_3_2.filter.rc.txt
