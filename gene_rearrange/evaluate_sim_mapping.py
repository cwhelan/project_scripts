#!/usr/local/bin/python

import sys
import pysam

true_locs_name = sys.argv[1]
mapped_locs_name = sys.argv[2]

true_locs = pysam.Samfile(true_locs_name,'rb')
mapped_locs = pysam.Samfile(mapped_locs_name,'rb')

lines_read = 0

while True:
    try:
        mapped_loc = mapped_locs.next()
        lines_read = lines_read + 1

        name = mapped_loc.qname
 #        print "read mapped " + name
        if mapped_loc.is_read1:
            mread1 = mapped_loc
        else:
            mread2 = mapped_loc
        mapped_loc = mapped_locs.next()
        lines_read = lines_read + 1
 #        print "read mapped name 2: " + mapped_loc.qname
        if mapped_loc.qname != name:
            print "assertion failed, second read name " + mapped_loc.qname + " is not equal to prior read " + name
        if mapped_loc.is_read1:
            mread1 = mapped_loc
        else:
            mread2 = mapped_loc

        if mread1.is_unmapped or mread2.is_unmapped:
            continue

        while True:
            try:
                true_loc = true_locs.next()
                tname = true_loc.qname
                if not tname == name:
                    true_loc = true_locs.next()
                    continue

                if true_loc.is_read1:
                    tread1 = true_loc
                else:
                    tread2 = true_loc
                true_loc = true_locs.next()
                if mapped_loc.qname != name:
                    print "assertion failed, second true read name " + mapped_loc.qname + " is not equal to prior read " + name

                if true_loc.is_read1:
                    tread1 = true_loc
                else:
                    tread2 = true_loc

                status = "WTF"
                if abs(abs(tread1.isize) - abs(mread1.isize)) <= 1:
                    status = "ISIZE_MATCH"
                else:
                    if abs(tread1.isize) > 500 and abs(mread1.isize) <= 500:
                        status = "FALSE_PE"
                    if abs(mread1.isize) > 500 and abs(tread1.isize) <= 500:
                        status = "FALSE_MP"
                    if abs(mread1.isize) > 500 and abs(tread1.isize) > 500:
                        status = "MISS_MP"
                    if abs(mread1.isize) < 500 and abs(tread1.isize) < 500:
                        status = "MISS_PE"            

                print "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}".format(name, tread1.pos, mread1.pos, tread2.pos, mread2.pos, abs(tread1.isize), abs(mread1.isize), status) 
                break
            except StopIteration:
                print "assertion failed: couldn't find the true location for read " + name
                sys.exit(-1)
    except StopIteration:
        break
print "Read {0} lines".format(lines_read)
        

