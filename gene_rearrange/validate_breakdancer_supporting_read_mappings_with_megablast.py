#!/usr/bin/env python

import os.path
import re
import subprocess
import sys
import tempfile

class alignment:
    def __init__(self, chrom, start, forward, record):
        self.chrom = chrom
        self.start = start
        self.forward = forward
        self.record = record

def concordant(loc1, loc2, mean_isize, mad_isize):
    if loc1.chrom != loc2.chrom:
        return False
    if loc1.start == loc2.start:
        return False
    if loc1.start < loc2.start and not (loc1.forward and not loc2.forward):
        return False
    if loc2.start < loc1.start and not (loc2.forward and not loc1.forward):
        return False
    if abs(loc1.start - loc2.start) > mean_isize + 10 * mad_isize:
        return False
    return True
    

def tr_chr_name(name):
    name = re.sub("^chr","",name)
    if (name == "M"):
        return "MT"
    if (name.find("gl") > -1):
        name = "GL" + name[name.find("gl") + 2:name.find("_random")] + ".1"
    return name

def lookup_sam_records(read_names, bam_file_name, region):
    read_name_tmp_file = tempfile.NamedTemporaryFile(mode='w', bufsize=0)
    read_name_tmp_file.write("\n".join(read_names) + "\n")

    samtools = subprocess.Popen(["/g/whelanch/software/bin/samtools", "view", bam_file_name, region], stdout=subprocess.PIPE)
    read_grep = subprocess.Popen(["grep", "-wf", read_name_tmp_file.name], stdin=samtools.stdout, stdout=subprocess.PIPE)
    samtools.stdout.close()
    sam_lines = read_grep.communicate()[0]
    read_name_tmp_file.close()
    return sam_lines

def sam_to_fasta(sam_lines):
    fasta_lines = ""
    for line in sam_lines.split("\n"):
        if line.strip() == "":
            continue
        fields = line.split()
        read_id = fields[0]
        seq = fields[9]
        flag = int(fields[1])
        if flag & 0x0040:
            read_num = "/1"
        else:
            read_num = "/2"
        fasta_line = ">" + read_id + read_num + "\n" + seq + "\n"
        fasta_lines += fasta_line
    return fasta_lines

def megablast(sam_lines, name, read_names, out_file, megablast_out_file):
    print "megablasting " + name
    fasta_tmp_file = tempfile.NamedTemporaryFile(mode='w', bufsize=0)
    fasta_tmp_file.write(sam_to_fasta(sam_lines))

    megablast_cmd = "megablast -d hg19 -i " + fasta_tmp_file.name + " -W 8 -G 6 -E 4 -F F -q -3 -r 2 -D 3 -s " + str(score_cutoff)
    megablast_output = subprocess.Popen(megablast_cmd.split(), stdout=subprocess.PIPE).communicate()[0]
    megablast_out_file.write("#\n# Executing MEGABLAST for supporting reads in " + name + "\n#\n")
    megablast_out_file.write(megablast_output)
    fasta_tmp_file.close()

    read1_alignments = {}
    read2_alignments = {}
    for line in megablast_output.split("\n"):
        if line.startswith("#") or line.strip() == "":
            continue
        fields = line.split()
        read_id = fields[0]
        [pair_id, pair_num] = read_id.split("/")
        chrom = fields[1]
        start = min(int(fields[8]), int(fields[9]))
        forward = int(fields[8]) < int(fields[9])
        loc = alignment(chrom, start, forward, line)
        if pair_num == "1":
            if not pair_id in read1_alignments:
                read1_alignments[pair_id] = []
            read1_alignments[pair_id].append(loc)
        else:
            if not pair_id in read2_alignments:
                read2_alignments[pair_id] = []
            read2_alignments[pair_id].append(loc)
    
    file_open = False
    read_names_with_concordant_mappings = set()
    for pair_id in read1_alignments.keys():
        r1_locs = read1_alignments[pair_id]
        if not pair_id in read2_alignments:
            continue
        r2_locs = read2_alignments[pair_id]
        found_concordant = False
        for r1_loc in r1_locs:
            if found_concordant:
                break
            for r2_loc in r2_locs:
                if concordant(r1_loc, r2_loc, mean_isize, mad_isize):
                    read_names_with_concordant_mappings.add(pair_id)
                    if not file_open:
                        concordant_mapping_file = open(name + "_concordant_mappings.txt", "w")
                        file_open = True
                    concordant_mapping_file.write("possible concordant mapping for " + pair_id + "\n")
                    concordant_mapping_file.write(r1_loc.record + "\n")
                    concordant_mapping_file.write(r2_loc.record + "\n")
                    found_concordant = True
                    break
    if file_open:
        concordant_mapping_file.close()
    print "\t".join([name, str(len(read_names)), str(len(read_names_with_concordant_mappings))])
    out_file.write("\t".join([name, str(len(read_names)), str(len(read_names_with_concordant_mappings))]) + "\n")

    fasta_tmp_file.close()


bd_read_file_name = sys.argv[1]
bam_file_name = sys.argv[2]
mean_isize = int(sys.argv[3])
mad_isize = int(sys.argv[4])
score_cutoff = int(sys.argv[5])
megablast_out_file_name = sys.argv[6]

bd_read_file_basename = os.path.basename(bd_read_file_name)

filter_chromosomes = False
if len(sys.argv) == 8:
    chromosome_filter = sys.argv[7] + "_"
    filter_chromosomes = True

name = ''
c1 = ''
s1 = 100000000000000
e1 = 0
c2 = ''
s2 = 100000000000000
e2 = 0

if not filter_chromosomes:
    out_file_name = bd_read_file_basename + ".validate.txt"
else:
    out_file_name = bd_read_file_basename + "_" + chromosome_filter + ".validate.txt"
out_file = open(out_file_name, "w")
bd_read_file = open(bd_read_file_name, "r")
megablast_out_file = open(megablast_out_file_name, "w")
for line in bd_read_file:
    if line.strip() == "":
        continue
    if line.startswith("track"):
        if name != '':
            if s2 < s1:
                ct = c1
                st = s1
                et = e1
                c1 = c2
                s1 = s2
                e1 = e2
                c2 = ct
                s2 = st
                e2 = et
            c1 = tr_chr_name(c1)
            c2 = tr_chr_name(c2)

            if not filter_chromosomes or name.startswith(chromosome_filter):
                # write out read names to a file
                # extract fasta from samtools view bam | grep -fw read name file
                r1 = c1 + ":" + str(s1) + "-" + str(e1)
                r2 = c2 + ":" + str(s2) + "-" + str(e2)
                #print "r1 = " + r1 + "; r2 = " + r2
                #print "\n".join(read_names) + "\n"

                sam_lines1 = lookup_sam_records(read_names, bam_file_name, r1)
                sam_lines2 = lookup_sam_records(read_names, bam_file_name, r2)

                sam_lines = sam_lines1 + sam_lines2
                
                megablast(sam_lines, name, read_names, out_file, megablast_out_file)

        name = ''
        c1 = ''
        s1 = 100000000000000
        e1 = 0
        c2 = ''
        s2 = 100000000000000
        e2 = 0
        read1 = True
        track_fields = line.split()
        name = track_fields[1].split("=")[1]
        read_names = set()
    else:
        line_fields = line.split()
        read_name = line_fields[3].split("|")[0]
        read_names.add(read_name)
        if read1:
            c1 = line_fields[0]
            if long(line_fields[1]) < s1:
                s1 = long(line_fields[1])
            if long(line_fields[2]) > e1:
                e1 = long(line_fields[2])
        else:
            c2 = line_fields[0]
            if long(line_fields[1]) < s2:
                s2 = long(line_fields[1])
            if long(line_fields[2]) > e2:
                e2 = long(line_fields[2])
        read1 = not read1


if s2 < s1:
    ct = c1
    st = s1
    et = e1
    c1 = c2
    s1 = s2
    e1 = e2
    c2 = ct
    s2 = st
    e2 = et
c1 = tr_chr_name(c1)
c2 = tr_chr_name(c2)

if not filter_chromosomes or name.startswith(chromosome_filter):
    # write out read names to a file
    # extract fasta from samtools view bam | grep -fw read name file
    r1 = c1 + ":" + str(s1) + "-" + str(e1)
    r2 = c2 + ":" + str(s2) + "-" + str(e2)
    #print "r1 = " + r1 + "; r2 = " + r2
    #print "\n".join(read_names) + "\n"

    sam_lines1 = lookup_sam_records(read_names, bam_file_name, r1)
    sam_lines2 = lookup_sam_records(read_names, bam_file_name, r2)
    sam_lines = sam_lines1 + sam_lines2
 
    megablast(sam_lines, name, read_names, out_file)

out_file.close()
megablast_out_file.close()
