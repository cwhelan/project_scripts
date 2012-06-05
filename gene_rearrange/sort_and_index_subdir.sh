#!/bin/bash

samtools sort -n novoalign.bam novoalign_sorted
samtools index novoalign_sorted.bam
