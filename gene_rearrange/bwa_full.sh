#!/bin/bash

/g/whelanch/software/bin/bwa aln -t 8 -n 2 /l2/users/whelanch/1000genomes/ref_bwa/human_g1k_v37.fasta $1 > $1.sai
/g/whelanch/software/bin/bwa aln -t 8 -n 2 /l2/users/whelanch/1000genomes/ref_bwa/human_g1k_v37.fasta $2 > $2.sai
/g/whelanch/software/bin/bwa sampe /l2/users/whelanch/1000genomes/ref_bwa/human_g1k_v37.fasta $1.sai $2.sai $1 $2 | /g/whelanch/software/bin/samtools view -Sb - > $3
