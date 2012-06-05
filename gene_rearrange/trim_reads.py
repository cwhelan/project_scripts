#!/usr/bin/env python                                                                                                                           

import sys
import gzip
from Bio import SeqIO

if sys.argv[1].endswith("gz"):
    in_handle = gzip.open(sys.argv[1], 'rb')
    out_handle = gzip.open(sys.argv[2], 'wb')
else:
    in_handle = open(sys.argv[1], 'r')
    out_handle = open(sys.argv[2], 'w')

trim_len = int(sys.argv[3])
trimmed_reads = (rec[:trim_len] for rec in \
                       SeqIO.parse(in_handle, "fastq"))
count = SeqIO.write(trimmed_reads, out_handle, "fastq")
print "Saved %i reads" % count

