#!/bin/bash
#samtools view -uF 2 /l2/users/whelanch/gene_rearrange/tier1_bwa.bam | \
#           bamToFastq -bam stdin \
#                      -fq1 /l2/users/whelanch/gene_rearrange/s.tier1.disc.1.fq \
#                      -fq2 /l2/users/whelanch/gene_rearrange/s.tier1.disc.2.fq
novoalign -d /l2/users/whelanch/1000genomes/ref_bwa/human_g1k_v37.nix -f /l2/users/whelanch/gene_rearrange/s_3_1_sequence.txt.tail
