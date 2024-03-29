#!/usr/bin/env python

import re
import subprocess
import os.path

script_path = os.path.dirname( os.path.realpath( __file__ ) )

prev_results_file = open("results.txt", "r")

first_line_fields = prev_results_file.readline().rstrip().split("\t")
bp_name = first_line_fields[2]
feature_name = first_line_fields[3]

num_real_bp_hits = first_line_fields[4]
prev_results_file.readline()
num_real_feature_hits = prev_results_file.readline().rstrip().split("\t")[4]
prev_results_file.readline()
real_bases_overlapped = prev_results_file.readline().rstrip().split("\t")[4]
prev_results_file.readline()

prev_results_file.close()

# plot the results
subprocess.call(map(str, ['Rscript', script_path + '/plot_results.R', feature_name, bp_name, num_real_bp_hits, num_real_feature_hits, real_bases_overlapped]))

