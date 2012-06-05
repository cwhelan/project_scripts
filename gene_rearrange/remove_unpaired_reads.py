#!/usr/local/bin/python

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


reads_discarded = 0
for read1_id, read1_seq, read1_qual in read1_iter:
    read2_id, read2_seq, read2_qual = read2_iter.next()
    
    r1_name = read1_id[0:read1_id.index('/')]
    r2_name = read2_id[0:read2_id.index('/')]

    while r1_name != r2_name:        
        reads_discarded += 1
        if (r1_name < r2_name):
            read1_id, read1_seq, read1_qual = read1_iter.next()
            r1_name = read1_id[0:read1_id.index('/')]
        else:
            read2_id, read2_seq, read2_qual = read2_iter.next()
            r2_name = read2_id[0:read2_id.index('/')]

    #SeqIO.write(read1, read1_out, "fastq")
    #SeqIO.write(read2, read2_out, "fastq")

    read1_out.write("@%s\n%s\n+\n%s\n" % (read1_id, read1_seq, read1_qual))
    read2_out.write("@%s\n%s\n+\n%s\n" % (read2_id, read2_seq, read2_qual))

print "Discarded {0} pairs\n".format(reads_discarded)

