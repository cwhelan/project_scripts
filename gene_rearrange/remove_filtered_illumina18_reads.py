#!/usr/bin/env python

# Usage python filter.py <read1.gz> <read2.gz>

# Takend from https://gist.github.com/1124940

from Bio import SeqIO
import itertools
import sys
import gzip
import os.path

f1 = gzip.open(sys.argv[1], "r")
f2 = gzip.open(sys.argv[2], "r")

head1, ext = os.path.splitext(sys.argv[1])
head2, ext = os.path.splitext(sys.argv[2])

f1o = gzip.open("%s_filtered.gz" % (head1), "wb")
f2o = gzip.open("%s_filtered.gz" % (head2), "wb")

itr1 = SeqIO.parse(f1, "fastq")
itr2 = SeqIO.parse(f2, "fastq")

for s1, s2 in itertools.izip(itr1, itr2):
    c1 = s1.description.split(":")
    c2 = s2.description.split(":")

    if c1[7] == "Y" or c2[7] == "Y":
        continue

    SeqIO.write([s1], f1o, "fastq")
    SeqIO.write([s2], f2o, "fastq")

f1.close()
f2.close()

f1o.close()
f2o.close()
