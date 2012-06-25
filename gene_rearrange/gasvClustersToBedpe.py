#!/usr/bin/env python

import sys

clusterFile = open(sys.argv[1], 'r')

for line in clusterFile:
    #print line
    if line.startswith("#"):
        #  print "starts with #"
        continue
    fields = line.split()
    cid = fields[0]
    c1 = fields[1]
    region1 = fields[2].split(",")
    s1 = region1[0]
    e1 = region1[0]
    c2 = fields[3]
    region2 = fields[4].split(",")
    #    print region2
    s2 = region2[0]
    e2 = region2[1]
    numPairs = fields[5]
    localization = fields[6]
    breaktype = fields[7]
    bedpeid = ";".join([cid, numPairs, localization, breaktype])
    o1 = "+"
    o2 = "-"
    if breaktype == "IR" or breaktype == "I+":
        o1 = "+"
        o2 = "+"
    elif breaktype == "I-":
        o1 = "-"
        o2 = "-"
    print "\t".join([c1, s1, e1, c2, s2, e2, bedpeid, numPairs, o1, o2])
    
