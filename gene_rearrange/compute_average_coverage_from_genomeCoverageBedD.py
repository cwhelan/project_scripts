#!/usr/local/bin/python

import sys
import math

chrom_length = 0
chrom_sum = 0

total_length = 0
total_sum = 0

read_any_lines = False
chrom = 'NA'

for line in sys.stdin:
    fields = line.strip().split()
    if read_any_lines and fields[0] != chrom:
        if chrom_length > 0:
            print "{0}\t{1}".format(chrom, chrom_sum / chrom_length)
        chrom_sum = 0
        chrom_length = 0
        
    read_any_lines = True
    chrom = fields[0]
    coverage = int(fields[3])
    length = 1
    if length == 0:
        continue
    coverage = float(coverage) * length
    
    total_sum += coverage
    chrom_sum += coverage

    chrom_length += length
    total_length += length

if (chrom_length > 0):
    print "{0}\t{1}".format(chrom, chrom_sum / chrom_length)
print "genome\t{0}".format(total_sum / total_length)    
