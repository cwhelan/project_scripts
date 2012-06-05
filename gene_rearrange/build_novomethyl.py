#!/usr/local/bin/python

import sys
import os
import shutil
import subprocess
import time

if not len(sys.argv) == 8:
	print "Usage: build_novomethyl.py working_dir reference read_file1 read_file2 target_isize isize_sd library_name"
	exit(1)

working_dir = sys.argv[1]
reference = sys.argv[2]
read_file1 = sys.argv[3]
read_file2 = sys.argv[4]
target_isize = sys.argv[5]
isize_sd = sys.argv[6]
library_name = sys.argv[7]

chunk_size = 500000

def split_file(name, prefix, chunk_size):
        split_cmd = "split -a 4 -d -l {2} {0} {1}.".format(name, prefix, chunk_size)
	print "executing: " + split_cmd
        subprocess.call(split_cmd, shell=True)

wcout = subprocess.Popen("wc -l {0}".format(read_file1), shell=True, stdout=subprocess.PIPE).communicate()[0]
file_length = int(wcout.split()[0])
print "Read file length: {0}".format(file_length)

num_chunks = (int(file_length) / chunk_size) + (1 if file_length % chunk_size > 0 else 0)
print "Number of chunks: {0}".format(num_chunks)

os.chdir(working_dir)
dagfile = open("{0}/novomethyl.dag".format(working_dir), 'w')
file1_base_name = os.path.basename(read_file1)
split_file(read_file1, file1_base_name, chunk_size)
file2_base_name = os.path.basename(read_file2)
split_file(read_file2, file2_base_name, chunk_size)

jobs = 0

first_alignment_chunk_job = jobs
for i in xrange(0,num_chunks):
	os.mkdir("wd{0}".format(i))
	shutil.move("{0}.{1:04d}".format(file1_base_name, i), "wd{0}".format(i))
	shutil.move("{0}.{1:04d}".format(file2_base_name, i), "wd{0}".format(i))
	chunk_job = i + jobs
	submit_file = "novomethyl.desc"
	dagfile.write("JOB {0} /l2/users/whelanch/gene_rearrange/submit_files/{2} DIR wd{1}/\n".format(chunk_job, i, submit_file))
	dagfile.write("VARS {0} read_file1=\"{2}/wd{1}/{3}.{1:04d}\"\n".format(chunk_job, i, working_dir, file1_base_name))
        dagfile.write("VARS {0} read_file2=\"{2}/wd{1}/{3}.{1:04d}\"\n".format(chunk_job, i, working_dir, file2_base_name))
	dagfile.write("VARS {0} reference=\"{1}\"\n".format(chunk_job, reference))
	dagfile.write("VARS {0} target_isize=\"{1}\"\n".format(chunk_job, target_isize))
	dagfile.write("VARS {0} isize_sd=\"{1}\"\n".format(chunk_job, isize_sd))
	dagfile.write("VARS {0} library_name=\"{1}\"\n".format(chunk_job, library_name))
	dagfile.write("RETRY {0} 2\n".format(chunk_job))

jobs = jobs + num_chunks

