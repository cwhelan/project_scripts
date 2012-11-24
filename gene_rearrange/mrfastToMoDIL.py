#!/usr/bin/env python

import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("insert_size", type=int, help="Insert size of the library")
args = parser.parse_args()

library_isize = args.insert_size

for line in sys.stdin:
    if line.startswith('@'):
        continue
    fields = line.split("\t")
    readid = fields[0]
    flag = int(fields[1])
    #only process first segments in template
    pos = int(fields[3])
    mpos = int(fields[7])
    if (mpos < pos):
        continue
    tlen = int(fields[8])
    sequence = fields[9]
    slen = len(sequence)
    endpos = pos + slen
    
    print "\t".join(map(str, [endpos, mpos, readid, library_isize - tlen, library_isize, readid]))
