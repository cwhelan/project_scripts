#!/usr/local/bin/python

import sys
import glob

samfile = open(sys.argv[1])
mbfiledir = sys.argv[2]
mbfile_glob = mbfiledir + "/*/*.mb"
mbfiles = glob.glob(mbfile_glob)

class ReadInfo:
	seq_name = ""
	seq = ""
	sam_mapping_chr = ""
	sam_mapping_loc = ""
	sam_concordant = False
	best_mb_mapping_chr = "--"
	best_mb_mapping_loc = "--"
	best_mb_escore = 0.0
	best_mb_bit_score = 0.0

read_infos = {}

while 1:
	line = samfile.readline().rstrip()
#	print "read line " + line
	if line == "":
		break
	fields = line.split()
	qname = fields[0]
#	print "read read " + qname
        seq = fields[9]
        flag = int(fields[1])
        if flag & 0x0040:
                name = qname + "/1"
        else:
                name = qname + "/2"
	read_info = ReadInfo()
	read_info.seq_name = name
	read_info.seq = seq
	read_info.sam_mapping_chr = fields[2]
	read_info.sam_mapping_loc = fields[3]
        if flag & 0x0002:
		read_info.sam_concordant = True
	read_infos[name] = read_info
#print "Read {0} sam records\n".format(len(read_infos.keys()))

for mbfilename in mbfiles:
#	print "Reading mb file " + mbfilename
	mbfile = open(mbfilename)
	while 1:	
		line = mbfile.readline().rstrip()
		if line == "":
			break
		if line.startswith('#'):
			continue
		fields = line.split()
		qname = fields[0]
		read_info = read_infos[qname]
		bitscore = float(fields[11])
		chrom = fields[1]
                loc1 = fields[8]
                loc2 = fields[9]
                left_loc = min(fields[8],fields[9])

		if (bitscore > read_info.best_mb_bit_score or (bitscore == read_info.best_mb_bit_score and chrom == read_info.sam_mapping_chr and left_loc == read_info.sam_mapping_loc)):
			escore = float(fields[10])
			read_info.best_mb_bit_score = bitscore
			read_info.best_mb_escore = escore
			read_info.best_mb_mapping_chr = chrom
			read_info.best_mb_mapping_loc = left_loc
	mbfile.close()

for qname in read_infos.keys():
	read_info = read_infos[qname]
	status = "MATCH"
	if (read_info.sam_mapping_chr != read_info.best_mb_mapping_chr) or (read_info.sam_mapping_loc != read_info.best_mb_mapping_loc):
		if read_info.sam_concordant:
			status = "CONCORDANT_CONFLICT"
		else:
			status = "DISCORDANT_CONFLICT"
	if read_info.best_mb_mapping_chr == "--":
		status = "NO_MEGABLAST_MAPPING"
	print "\t".join([read_info.seq_name, read_info.seq, read_info.sam_mapping_chr, read_info.sam_mapping_loc, read_info.best_mb_mapping_chr, read_info.best_mb_mapping_loc, str(read_info.best_mb_escore),str(read_info.best_mb_bit_score), status])

		
