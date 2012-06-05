#!/usr/local/bin/python

import sys

total_depth = 0
len = 0

for line in sys.stdin:
	fields = line.split()
	depth = long(fields[3])
	total_depth += depth
	len = len + 1

print "Total base reads: {0}".format(total_depth)
print "Average coverage: {0}".format((total_depth / len))
