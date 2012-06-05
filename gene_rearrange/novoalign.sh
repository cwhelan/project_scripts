#!/bin/bash
#samtools view -uF 2 /l2/users/whelanch/gene_rearrange/tier1_bwa.bam | \
#           bamToFastq -bam stdin \
#                      -fq1 /l2/users/whelanch/gene_rearrange/s.tier1.disc.1.fq \
#                      -fq2 /l2/users/whelanch/gene_rearrange/s.tier1.disc.2.fq
/g/whelanch/software/bin/novoalign -d /l2/users/whelanch/1000genomes/ref_bwa/human_g1k_v37.nix \
                         -f /l2/users/whelanch/gene_rearrange/s_3_1_sequence.txt /l2/users/whelanch/gene_rearrange/s_3_2_sequence.txt \
                         -F ILMFQ -i MP 3000-8000 250,100 -a -r None -R 0 -oSAM $'@RG\tID:readgroup\tPU:platform­ unit\tLB:library' | \
	/g/whelanch/software/bin/samtools view -Sb - > /l2/users/whelanch/gene_rearrange/novoalign.bam 

