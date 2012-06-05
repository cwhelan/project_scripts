#!/usr/local/bin/python

import sys

for line in sys.stdin:
	fields = line.split()
	mate_chr = fields[6]
	flag = int(fields[1])
	if (not flag & 0x2):
		continue
	if (mate_chr != '='):
		continue
	chr = fields[2]
	pos = int(fields[3])
	mate_pos = int(fields[7])

	isize = int(fields[8])

	if (flag & 0x001):
		

	print "\t".join([chr, str(start), str(end), str(end - start), fields[8]])
	
