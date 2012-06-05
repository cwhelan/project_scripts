#!/bin/bash

#samtools view -F 0x000E -b /l2/users/whelanch/gene_rearrange/novoalign.bam > /l2/users/whelanch/gene_rearrange/novo_discordants.bam
#bamToBed -i /l2/users/whelanch/gene_rearrange/novo_discordants.bam -tag NM > /l2/users/whelanch/gene_rearrange/novo_discordants.bed
#pairDiscordants.py -i /l2/users/whelanch/gene_rearrange/novo_discordants.bed -m hydra -y 3000 -z8000 > /l2/users/whelanch/gene_rearrange/novo_discordants.bedpe
#dedupDiscordants.py -i /l2/users/whelanch/gene_rearrange/novo_discordants.bedpe -s 3 > /l2/users/whelanch/gene_rearrange/novo_discordants.dedup.bedpe
#hydra -in /l2/users/whelanch/gene_rearrange/novo_discordants.dedup.bedpe -out /l2/users/whelanch/gene_rearrange/novoalign.breaks -mld 8120 -mno 21165 
hydra -in /l2/users/whelanch/gene_rearrange/novo_discordants.dedup.bedpe -out /l2/users/whelanch/gene_rearrange/novoalign_ms3.breaks -mld 8120 -mno 21165 -ms 3
