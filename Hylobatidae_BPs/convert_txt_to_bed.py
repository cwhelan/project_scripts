#!/usr/bin/env python

import sys

for line in sys.stdin:
    line = line.replace(":", "\t")
    line = line.replace("-", "\t")
    line = line.replace(",", "")
    fields = line.rstrip().split("\t")
    start = min(int(fields[1]), int(fields[2]))
    end = max(int(fields[1]), int(fields[2]))
    print "\t".join(map(str, [fields[0], start, end, fields[3].replace(" ", ""), "1"]))
