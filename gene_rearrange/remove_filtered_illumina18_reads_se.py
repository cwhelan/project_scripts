#!/usr/bin/env python

# Usage python filter.py <read1.gz> <read2.gz>

# Takend from https://gist.github.com/1124940

from Bio import SeqIO
import itertools
import sys
import gzip
import os.path

f1 = gzip.open(sys.argv[1], "r")

head1, ext = os.path.splitext(sys.argv[1])

f1o = gzip.open("%s_filtered.gz" % (head1), "wb")

itr1 = SeqIO.parse(f1, "fastq")

reads_filtered = 0
for s1 in itr1:
    c1 = s1.description.split(":")

    if c1[7] == "Y":
        reads_filtered = reads_filtered + 1
        continue

    SeqIO.write([s1], f1o, "fastq")

f1.close()

f1o.close()

print "removed " + str(reads_filtered) + " reads"
