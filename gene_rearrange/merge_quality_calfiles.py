#!/usr/bin/env python

import sys
import glob

working_dir = sys.argv[1]

file_pattern = working_dir + "/wd*/calfile.txt"

calfiles = glob.glob(file_pattern)

qualities = {}

for calfilename in calfiles:
    calfile = open(calfilename)
    header = True
    for line in calfile:
        if line.startswith("#"):
            continue
        if header:
            header = False
            continue
        fields = line.strip().split()
        key = tuple([int(fields[0]), int(fields[1]), int(fields[2]), fields[3]])
        val = tuple(map(int, fields[4:len(fields)]))
        if not key in qualities:
            qualities[key] = val
        else:
            qualities[key] = map(sum, zip(val, qualities[key]))

output = open(working_dir + "/calfile.txt", "w")

output.write("Side\tIndex\tQuality\tBase\tA's\tC's\tG's\tT's\tNMs\n")
keys = qualities.keys()
keys.sort()
for key in keys:
    output.write("\t".join(map(str, list(key) + list(qualities[key]))))
    output.write("\n")
