#!/usr/local/bin/python

import sys
import pysam

bam_1 = sys.argv[1]
bam_2 = sys.argv[2]

s1 = pysam.Samfile(bam_1,'rb')
s2 = pysam.Samfile(bam_2,'rb')

pairs = 0
inter_chrom = 0
short_inserts = 0
long_inserts = 0
sum_long_isize = 0
sum_short_isize = 0

while True:
    try:
        r1 = s1.next()
        r2 = s2.next()
        if not r1.is_unmapped and not r2.is_unmapped:
            pairs += 1
            if not r1.tid == r2.tid:
                inter_chrom += 1
                continue
            long_insert = (abs(r1.pos - r2.pos) > 500)
            if not long_insert:
                continue
            if long_insert:
                long_inserts += 1
            else:
                short_inserts += 1
            if r2.is_reverse and long_insert:
                leftread = r2
                rightread = r1
            else:
                leftread = r1
                rightread = r2
            chrom = s1.getrname(r1.tid)
            lpos = leftread.pos
            rpos = rightread.pos + rightread.rlen
            isize = rpos - lpos
            if long_insert:
                sum_long_isize += isize
            else:
                sum_short_isize += isize
            print "{0}\t{1}\t{2}\t{3}".format(chrom, lpos, rpos, isize)
            
            
    except StopIteration:
        break
    
sys.stderr.write("read {0} pairs\n".format(pairs))
sys.stderr.write("inter chromosome: {0}\n".format(inter_chrom))
sys.stderr.write("mean long isize: {0}\n".format(sum_long_isize/(long_inserts - inter_chrom)))
sys.stderr.write("mean short isize: {0}\n".format(sum_short_isize/(short_inserts - inter_chrom)))

