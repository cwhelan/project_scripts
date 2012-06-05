#!/usr/bin/python

import sys

for line in sys.stdin:
    if line[0] == '#':
        continue
    fields = line.rstrip().split()
    chrom = fields[0]
    start = int(fields[1]) - 1
    tags = {}
#    print "tags: " + fields[7]
    for tag in fields[7].split(";"):
        tag_fields = tag.split("=")
        k = tag_fields[0]
        if len(tag_fields) == 2:
            v = tag_fields[1]
        else:
            v = True
        tags[k] = v
    var_id = "NA"
    if ("NOVEL" in tags):
        var_id = "NOVEL"
    if ("DBVARID" in tags):
        var_id = tags["DBVARID"]
    svlen = -1 * int(tags["SVLEN"])
    end = start + svlen
    print "\t".join([chrom, str(start), str(end), var_id, str(svlen)])
    
    
