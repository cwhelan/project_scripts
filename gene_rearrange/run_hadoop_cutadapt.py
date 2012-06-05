#!/usr/bin/python

import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.QualityIO import FastqGeneralIterator

#read1_iter = SeqIO.parse(sys.argv[1], "fastq")
#read2_iter = SeqIO.parse(sys.argv[2], "fastq")

read1_iter = FastqGeneralIterator(open(sys.argv[1]))
read2_iter = FastqGeneralIterator(open(sys.argv[2]))


kv_out_file = open(sys.argv[1] + ".tmp", "w")

for read1 in read1_iter:
    read2 = read2_iter.next()
    
    print (read1[0])

    # strip off the /1 in read1 - this will be the key for the MR data file
    read_id = read1[0][:len(read1[0])-2]
    
    kv_out_file.write("\t".join([read_id, read1[0], read1[1], read1[2], read2[0], read2[1], read2[2]]) + "\n")
    



