#!/usr/local/bin/python

import sys

for line in sys.stdin:
	if (line.startswith("@")):
		print line.rstrip()
		continue
	fields = line.split()
	if fields[8] == '*':
		print line.rstrip()
	insert_size = int(fields[8])
	if (insert_size == 0 or abs(insert_size) > 500):
		print line.rstrip()
