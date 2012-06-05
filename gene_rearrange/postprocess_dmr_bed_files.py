#!/usr/bin/env python

import sys

track_name = sys.argv[1]
hypo_regions = open(sys.argv[2], 'r')
hyper_regions = open(sys.argv[3], 'r')
insignificant_regions = open(sys.argv[4], 'r')

print "track name=\"" + track_name + "\" itemRgb=\"On\""

for line in hypo_regions:
    fields = line.rstrip().split()
    print "\t".join([fields[0], fields[1], fields[2], fields[3], fields[4], fields[5], fields[1], fields[2], "0,0,255"])

for line in hyper_regions:
    fields = line.rstrip().split()
    print "\t".join([fields[0], fields[1], fields[2], fields[3], fields[4], fields[5], fields[1], fields[2], "255,0,0"])

for line in insignificant_regions:
    fields = line.rstrip().split()
    print "\t".join([fields[0], fields[1], fields[2], fields[3], fields[4], fields[5], fields[1], fields[2], "190,190,190"])

