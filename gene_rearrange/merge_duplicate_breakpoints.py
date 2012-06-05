#!/usr/bin/env python

import sys

class BEDPE_HYDRA (object):
	"""
	BEDPE class.  Can initialize to "NULL"
	by passign an empty list
	"""	
	def __init__(self, bedList = []):
		if len(bedList) > 0:
			self.lineNum   = bedList[0]
			self.c1        = bedList[1]
			self.s1        = int(bedList[2])
			self.e1        = int(bedList[3])
			self.c2        = bedList[4]
			self.s2        = int(bedList[5])
			self.e2        = int(bedList[6])
			self.name      = bedList[7]
			self.score     = bedList[8]
			self.o1        = bedList[9]		
			self.o2        = bedList[10]
			self.edit1     = float(bedList[11])		
			self.edit2     = float(bedList[12])

			self.totalEdit     = self.edit1 + self.edit2
			self.valid = 1
		else:
			self.valid = 0
			
	def writeBed(self):
		print self.c1 + "\t" + str(self.s1) + "\t" + str(self.e1) + "\t" + \
			self.c2 + "\t" + str(self.s2) + "\t" + str(self.e2) + "\t" + \
			self.name + "\t" + self.score + "\t" + \
			self.o1 + "\t" + self.o2 + "\t" + \
			str(self.edit1) + "\t" + str(self.edit2)

        def shortStr(self):
            return self.c1 + ":" + str(self.s1) + "-" + str(self.e1) + "/" + self.c2 + ":" + str(self.s2) + "-" + str(self.e2)

bedpeIn = sys.argv[1]
slop = int(sys.argv[2])

all_lines = []
line_num = 0
for line in open(bedpeIn, 'r'):
    line = str(line_num) + "\t" + line
    line_num = line_num + 1
    lineList = line.strip().split()
    if (len(lineList) > 0):
        # create a BEDPE from the current line and if there
        # is a valid previous line, then check for duplicates
        all_lines.append(BEDPE_HYDRA(lineList))
        

line_num = 0
# for line in open(bedpeIn, 'r'):
#     line = str(line_num) + "\t" + line
#     line_num = line_num + 1
#     lineList = line.strip().split()
#     if (len(lineList) > 0):
#         # create a BEDPE from the current line and if there
#         # is a valid previous line, then check for duplicates
for curr in all_lines:
    for bpe in all_lines:
        if curr.lineNum == bpe.lineNum:
            continue
#        if bpe.c1 > curr.c1:
#            break
        if curr.c1 != bpe.c1:
            continue
        if (abs(curr.s1 - bpe.e1) < slop or abs(curr.e1 - bpe.s1) < slop):
            if (curr.c2 == bpe.c2) and (abs(curr.s2 - bpe.e2) < slop or abs(curr.e2 - bpe.s2) < slop):
                print "removing " + bpe.shortStr() + " as a dupe of " + curr.shortStr()
                all_lines.remove(bpe)

for bpe in all_lines:
    bpe.writeBed()
