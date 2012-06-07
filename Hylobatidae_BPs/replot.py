#!/usr/bin/env python

import re

prev_results_file = open("results.txt", "r")

first_line_fields = prev_results_file.readline().rstrip().split("\t")
description = first_line_fields[0]
description = description.replace("_", " ")
bp_name = re.match("Number of (.+)s that Overlap a (.+)$", s).groups()[0]
feature_name = re.match("Number of (.+)s that Overlap a (.+)$", s).groups()[1]

num_real_bp_hits = first_line_fields[2]
prev_results_file.readline()
num_real_feature_hits = prev_results_file.readline().rstrip().split("\t")[2]
prev_results_file.readline()
num_real_bases_overlapped = prev_results_file.readline().rstrip().split("\t")[2]
prev_results_file.readline()

prev_results_file.close()

# plot the results
subprocess.call(map(str, ['/g/whelanch/software/bin/Rscript', script_path + '/plot_results.R', feature_name, bp_name, num_real_bp_hits, num_real_feature_hits, real_bases_overlapped]))

