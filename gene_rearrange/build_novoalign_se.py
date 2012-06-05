#!/usr/local/bin/python

import sys
import os
import shutil
import subprocess
import time

if not len(sys.argv) == 6:
	print "Usage: build_novoalign_se.py working_dir reference read_file repeat_report library_name" 
	exit(1)

working_dir = sys.argv[1]
reference = sys.argv[2]
read_file = sys.argv[3]
repeat_report = sys.argv[4]
library_name=sys.argv[5]

chunk_size = 1000000

def split_file(name, prefix, chunk_size):
        split_cmd = "split -a 4 -d -l {2} {0} {1}.".format(name, prefix, chunk_size)
        subprocess.call(split_cmd, shell=True)

wcout = subprocess.Popen("wc -l {0}".format(read_file), shell=True, stdout=subprocess.PIPE).communicate()[0]
file_length = int(wcout.split()[0])
num_chunks = (int(file_length) / chunk_size) + 1

os.chdir(working_dir)
dagfile = open("{0}/novoalign_first_tier.dag".format(working_dir), 'w')
file_base_name = os.path.basename(read_file).split(".")[0]
split_file(read_file, file_base_name, chunk_size)

jobs = 0

first_alignment_chunk_job = jobs
for i in xrange(jobs,num_chunks+jobs):
	os.mkdir("wd{0}".format(i))
	shutil.move("{0}.{1:04d}".format(file_base_name, i), "wd{0}".format(i))
	dagfile.write("JOB {0} /l2/users/whelanch/gene_rearrange/submit_files/novoalign_se.desc DIR wd{0}/\n".format(i))
	dagfile.write("VARS {0} read_file=\"{1}/wd{0}/{2}.{0:04d}\"\n".format(i, working_dir, file_base_name))
	dagfile.write("VARS {0} reference=\"{1}\"\n".format(i, reference))
	dagfile.write("VARS {0} repeat_report=\"{1}\"\n".format(i, repeat_report))
	dagfile.write("VARS {0} library_name=\"{1}\"\n".format(i, library_name))

jobs = jobs + num_chunks

merge_job = jobs
dagfile.write("JOB {0} /l2/users/whelanch/gene_rearrange/submit_files/merge_bams.desc DIR {1}\n".format(merge_job, working_dir))
dagfile.write("VARS {0} output=\"novoalign_tier1_nsort\"\n".format(merge_job))
dagfile.write("VARS {0} working_dir=\"{1}\"\n".format(merge_job, working_dir))
dagfile.write("VARS {0} in_file_names=\"{1}_sorted.bam\"\n".format(merge_job, file_base_name))

jobs = jobs + 1

cleanup_job = jobs
dagfile.write("JOB {0} /l2/users/whelanch/gene_rearrange/submit_files/cleanup_working_dirs.desc DIR {1}\n".format(cleanup_job, working_dir))
dagfile.write("VARS {0} working_dir=\"{1}\"\n".format(cleanup_job, working_dir))


dagfile.write("PARENT ")
for i in xrange(first_alignment_chunk_job,num_chunks+first_alignment_chunk_job):
	dagfile.write("{0} ".format(i))
dagfile.write("CHILD {0}\n".format(num_chunks+first_alignment_chunk_job))

dagfile.write("PARENT {0} CHILD {1}\n".format(merge_job, cleanup_job))
