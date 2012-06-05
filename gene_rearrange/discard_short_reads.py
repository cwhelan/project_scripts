#!/usr/bin/env python

# Throws away pairs of reads in which one of the reads is less than
# min_length, and trims both pairs of reads to max_length.
#
# Author: Chris Whelan

import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.QualityIO import FastqGeneralIterator

#read1_iter = SeqIO.parse(sys.argv[1], "fastq")
#read2_iter = SeqIO.parse(sys.argv[2], "fastq")

read1_iter = FastqGeneralIterator(open(sys.argv[1]))
read2_iter = FastqGeneralIterator(open(sys.argv[2]))

read1_out = open(sys.argv[3], 'w')
read2_out = open(sys.argv[4], 'w')

min_length = int(sys.argv[5])
max_length = int(sys.argv[6])

pairs_discarded = 0
for read1_id, read1_seq, read1_qual in read1_iter:
    read2_id, read2_seq, read2_qual = read2_iter.next()
    
    read1_seq = read1_seq[5:]
    read1_qual = read1_qual[5:]
    read2_seq = read2_seq[5:]
    read2_qual = read2_qual[5:]
    
    if len(read1_seq) < min_length or len(read2_seq) < min_length:
        pairs_discarded += 1
        continue
    #SeqIO.write(read1, read1_out, "fastq")
    #SeqIO.write(read2, read2_out, "fastq")

    if len(read1_seq) > max_length:
        read1_seq = read1_seq[:max_length]
        read1_qual = read1_qual[:max_length]

    if len(read2_seq) > max_length:
        read2_seq = read2_seq[:max_length]
        read2_qual = read2_qual[:max_length]

    read1_out.write("@%s\n%s\n+\n%s\n" % (read1_id, read1_seq, read1_qual))
    read2_out.write("@%s\n%s\n+\n%s\n" % (read2_id, read2_seq, read2_qual))

print "Discarded {0} pairs\n".format(pairs_discarded)

