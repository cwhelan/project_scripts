#!/usr/bin/python

import sys
import random
import subprocess

wcout = subprocess.Popen("samtools view -F 4 " + sys.argv[1] + " | wc -l ", shell=True, stdout=subprocess.PIPE).communicate()[0]
lines = int(wcout.split()[0])
selections = frozenset(random.sample(xrange(lines), int(sys.argv[2])))
max_sel = max(selections)
#print selections
#sys.stderr.write("lines: %d\n" % lines)
f = subprocess.Popen("samtools view -F 4 " + sys.argv[1], shell=True, stdout=subprocess.PIPE)
l = 0
while 1:
	l1 = f.stdout.readline()
	if l1 == "":
		break
	if i > max_sel:
		break
	if l in selections:
		print l1.rstrip()
	l = l + 1
