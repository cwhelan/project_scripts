#!/usr/local/bin/python

import sys

for line in sys.stdin:
	if line.startswith('@'):
		continue
	fields = line.strip().split()
	if fields[9] == '*':
		print fields[0]

