#!/bin/bash

samtools view -uF 2 novoalign.bam | \
           bamToFastq -bam stdin \
           -fq1 s_3_1.novo1.disc.1.fq \
           -fq2 s_3_2.novo1.disc.2.fq
novoalign -d /l2/users/whelanch/1000genomes/ref_bwa/human_g1k_v37.nix \
                         -f s_3_1.novo1.disc.1.fq s_3_2.novo1.disc.2.fq \
                         -i MP 3000-8000 250,100 -r Ex 500 -t 180 -o SAM > novoalign_sensitive.sam
cat novoalign_sensitive.sam | samtools view -Sb - > novoalign_sensitive.bam 

