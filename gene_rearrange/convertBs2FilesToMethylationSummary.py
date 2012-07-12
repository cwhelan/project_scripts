#!/usr/bin/env python

# convert's Larry's .bs2 files into the _methylation_summary.txt.gz needed by my R scripts

import sys
import glob
import os
import gzip

wd = sys.argv[1]
prefix = sys.argv[2]

outfile = gzip.open(prefix + "_methylation_summary.txt.gz", "wb")

# bs2 fields:
#<POS>   <COV_DEPTH>     <ME_RATE>       <TOP_RATE>      <BOT_RATE>      <TOP_M> <TOP_UM>        <BOT_M> <BOT_UM>
#77942   36      0.61    0.60    0.61    3       2       19      12

for fname in glob.glob(wd + "/*.bs2"):
    print "working on " + os.path.basename(fname)
    chrom = os.path.basename(fname).split(".")[0]
    f = open(fname, "r")
    f.readline()
    for line in f.readlines():
        fields = line.rstrip().split()
        pos = fields[0]
        meth = int(fields[5]) + int(fields[7])
        unmeth = int(fields[6]) + int(fields[8])
        assert meth + unmeth == int(fields[1])
        outfile.write("\t".join([chrom, str(pos), "+", str(meth), str(unmeth)]) + "\n")
    f.close()
outfile.close()
                    
