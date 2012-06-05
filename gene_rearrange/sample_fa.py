#!/usr/bin/python

import sys
import random
import subprocess

wcout = subprocess.Popen("wc -l " + sys.argv[1], shell=True, stdout=subprocess.PIPE).communicate()[0]
lines = int(wcout.split()[0]) / 2
selections = frozenset(random.sample(xrange(lines), int(sys.argv[2])))
#print selections
#sys.stderr.write("lines: %d\n" % lines)
f = open(sys.argv[1])
l = 0
while 1:
	#if l % 100000 == 0:
	#	sys.stderr.write("%d\n" % l)
	l1 = f.readline()
	#print "read l1 = " + l1
	if l1 == "":
		break
	l2 = f.readline()
	#print "read l2 = " + l2
	#print l
	if l in selections:
		print l1.rstrip()
		print l2.rstrip()
	l = l + 1
