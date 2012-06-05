#!/usr/bin/env python

# Author: Christopher Whelan (whelanch@ohsu.edu)
#
# This script calculates the methylation call rate at each read position to help identify biases
# in sequencing from Bismark output.

import sys
import gzip
import re

from numpy import *

vanilla = False
if len(sys.argv) != 4:
    if len(sys.argv) == 5 and sys.argv[4] == "--vanilla":
        vanilla = True
    else:
        print "Usage: calculateBismarkBias.py bismark_pe_mapping_report max_read_length library_name [--vanilla]"
        sys.exit()

if sys.argv[1].endswith("gz"):
    bismark_output = gzip.open(sys.argv[1], "r")
else:
    bismark_output = open(sys.argv[1], "r")
max_read_length = int(sys.argv[2])
library_name = sys.argv[3]

r1_methylated_calls_by_position = array([0] * max_read_length, dtype=float)
r1_total_calls_by_position = array([0] * max_read_length, dtype=float)

r2_methylated_calls_by_position = array([0] * max_read_length, dtype=float)
r2_total_calls_by_position = array([0] * max_read_length, dtype=float)

def increment_counts(call_string, methylated_calls_by_position, total_calls_by_position):
    pos = 0
    for c in call_string:
        if c != '.':
            if c.isupper():
                methylated_calls_by_position[pos] += 1
            total_calls_by_position[pos] += 1
        pos += 1

sam_pattern = re.compile("XM:Z:(\S+)")
num_lines = 0
for line in bismark_output:
    if vanilla:
        if line.startswith("Bismark"):
            continue
        num_lines += 1
        fields = line.strip().split()
        call_string_1 = fields[8]
        call_string_2 = fields[11]
        increment_counts(call_string_1, r1_methylated_calls_by_position, r1_total_calls_by_position)
        increment_counts(call_string_2, r2_methylated_calls_by_position, r2_total_calls_by_position)
    else:
        if line.startswith("@"):
            continue
        num_lines += 1
        call_string = sam_pattern.search(line).group(1)
        flag = int(line.split("\t")[1])
        if flag & 0x0040:
            increment_counts(call_string, r1_methylated_calls_by_position, r1_total_calls_by_position)
        elif flag & 0x0080:
            increment_counts(call_string, r2_methylated_calls_by_position, r2_total_calls_by_position)
        else:
            print "Bad line!!! " + line
        

r1_pcts = r1_methylated_calls_by_position / r1_total_calls_by_position
r1_pcts = nan_to_num(r1_pcts)

r2_pcts = r2_methylated_calls_by_position / r2_total_calls_by_position
r2_pcts = nan_to_num(r2_pcts)

#print r1_methylated_calls_by_position
#print r1_total_calls_by_position
print "Examined " + str(num_lines) + " alignments"
print "Read 1 methylation calls by position:"
print "\t".join(map(str, r1_pcts))
print "Read 2 methylation calls by position:"
print "\t".join(map(str, r2_pcts))

R_plot_file = open(library_name + "_bias_by_position.R", "w")

R_plot_file.write("r1 <- c(" + ", ".join(map(str, r1_pcts)) + ")\n")
R_plot_file.write("r2 <- c(" + ", ".join(map(str, r2_pcts)) + ")\n")
R_plot_file.write("pdf(\"" + library_name + "_bias_by_position.pdf\")\n")
R_plot_file.write("par(mfrow=c(2,1), oma = c( 0, 0, 0, 0 ))\n")
R_plot_file.write("plot(seq(r1), 100 * r1, type=\"o\", xlab=\"Position in Read\", ylab=\"% C called methylated\", ylim=c(0,7), xlim=c(0,"+str(max_read_length)+"), main=\"Read 1\")\n")
R_plot_file.write("plot(seq(r2), 100 * r2, type=\"o\", xlab=\"Position in Read\", ylab=\"% C called methylated\", ylim=c(0,7), xlim=c(0,"+str(max_read_length)+"), main=\"Read 2\")\n")
R_plot_file.write("dev.off()\n")

R_plot_file.close()
print "plot the results by executing \"Rscript " + library_name + "_bias_by_position.R\""
    
