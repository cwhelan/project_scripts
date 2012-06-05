#!/usr/local/bin/python

import sys
import os
import shutil
import subprocess
import glob
import time

working_dir = sys.argv[1]
reference = sys.argv[2]
read_file1 = sys.argv[3]
read_file2 = sys.argv[4]
ignore_suffix = sys.argv[5]

chunk_size = 800000

def split_file(name, prefix, chunk_size):
	split_cmd = "split -a 4 -d -l {2} {0} {1}.".format(name, prefix, chunk_size)
	subprocess.Popen(split_cmd, shell=True)

wcout = subprocess.Popen("wc -l {0}".format(read_file1), shell=True, stdout=subprocess.PIPE).communicate()[0]
file_length = int(wcout.split()[0])
num_chunks = (int(file_length) / chunk_size) + 1

print "Processing file {0} into {1} chunks of size {2}\n".format(read_file1, num_chunks, chunk_size)
	
os.chdir(working_dir)
dagfile = open("{0}/pash.dag".format(working_dir), 'w')
file1_base_name = os.path.basename(read_file1)
split_file(read_file1, file1_base_name, chunk_size)
file2_base_name = os.path.basename(read_file2)
split_file(read_file2, file2_base_name, chunk_size)

# sleeping a bit to let the file system catch up
time.sleep(60)

for i in xrange(0,num_chunks):
	os.mkdir("{0}".format(i))
	shutil.move("{2}/{0}.{1:04d}".format(file1_base_name, i, working_dir), "{0}/{1}".format(working_dir, i))
	shutil.move("{2}/{0}.{1:04d}".format(file2_base_name, i, working_dir), "{0}/{1}".format(working_dir, i))
	dagfile.write("JOB {0} /l2/users/whelanch/gene_rearrange/submit_files/pash_dag.desc DIR {1}/{2}/\n".format(i*2, sys.argv[1], i))
	dagfile.write("VARS {0} inputfile=\"{1}.{2:04d}\"\n".format(i*2, file1_base_name, i))
	dagfile.write("VARS {0} reference=\"{1}\"\n".format(i*2, reference))
	dagfile.write("VARS {0} ignoresuffix=\"{1}\"\n".format(i*2, ignore_suffix))
	dagfile.write("JOB {0} /l2/users/whelanch/gene_rearrange/submit_files/pash_dag.desc DIR {1}/{2}/\n".format(i*2 + 1, sys.argv[1], i))
	dagfile.write("VARS {0} inputfile=\"{1}.{2:04d}\"\n".format(i*2 + 1, file2_base_name, i))
	dagfile.write("VARS {0} reference=\"{1}\"\n".format(i*2 + 1, reference))
	dagfile.write("VARS {0} ignoresuffix=\"{1}\"\n".format(i*2 + 1, ignore_suffix))

