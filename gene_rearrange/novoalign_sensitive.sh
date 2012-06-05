#!/bin/bash

samtools view -uF 2 /l2/users/whelanch/gene_rearrange/novoalign.bam | \
           bamToFastq -bam stdin \
           -fq1 /l2/users/whelanch/gene_rearrange/s_3_1.novo1.disc.1.fq \
           -fq2 /l2/users/whelanch/gene_rearrange/s_3_2.novo1.disc.2.fq
novoalign -c 4 -d /l2/users/whelanch/1000genomes/ref_bwa/human_g1k_v37.nix -a \
                         -f /l2/users/whelanch/gene_rearrange/s_3_1.novo1.disc.1.fq /l2/users/whelanch/gene_rearrange/s_3_2.novo1.disc.2.fq \
                         -i 500 50 -r Ex 1100 -t 300 -o SAM | \
            samtools view -Sb - > /l2/users/whelanch/gene_rearrange/novoalign_tier2.bam 

