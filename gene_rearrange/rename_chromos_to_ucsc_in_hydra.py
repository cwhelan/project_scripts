#!/usr/bin/python

import sys

def tr_name(name):
	if (name == "MT"):
		return "chrM"
	if (name.startswith("GL")):
		gl_num = int(name[2:name.find('.')])
		if (gl_num >= 193 and gl_num <= 194):
			return "chr4_gl%(num)06d_random" % {"num": gl_num}	
		if (gl_num == 195):
			return "chr7_gl000195_random"
		if (gl_num >= 203 and gl_num <= 206):
			return "chr17_gl%(num)06d_random" % {"num": gl_num}
		if (gl_num == 207):
			return "chr18_gl%(num)06d_random" % {"num": gl_num}
                if (gl_num >= 208 and gl_num <= 209):
                        return "chr19_gl%(num)06d_random" % {"num": gl_num}
		if (gl_num >= 211):
			return "chrUn_gl%(num)06d" % {"num": gl_num}
	return "chr" + name

for line in sys.stdin:
	fields = line.rstrip().split()
	fields[0] = tr_name(fields[0])
	fields[3] = tr_name(fields[3])
	print "\t".join(fields)
