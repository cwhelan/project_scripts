#!/usr/bin/env python

import sys

meth = 0
unmeth = 0

curr_chrom = "NA"
curr_loc = "0"
curr_strand = "NA"

for line in sys.stdin.readlines():
    if line.startswith("Bismark"):
        continue
    fields = line.rstrip().split("\t")
    if len(fields) < 5:
        sys.stderr.write('bad line' + line + '\n')
    chrom = fields[2]
    loc = fields[3]
    strand = fields[1]
    call = fields[4]
    if (chrom != curr_chrom or loc != curr_loc):
        if curr_chrom != "NA":
            print "\t".join([curr_chrom, curr_loc, curr_strand, str(meth), str(unmeth)])
        meth = 0
        unmeth = 0
        curr_chrom = chrom
        curr_loc = loc
        curr_strand = strand        
    if (call.isupper()):
        meth += 1
    else:
        unmeth += 1
print "\t".join([curr_chrom, curr_loc, curr_strand, str(meth), str(unmeth)])
        
