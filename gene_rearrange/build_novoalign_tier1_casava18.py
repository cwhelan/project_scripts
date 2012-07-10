#!/usr/bin/env python

import sys
import os
import shutil
import subprocess
import time
import glob

if not len(sys.argv) == 13:
	print "Usage: build_novoalign_tier1.py working_dir reference read_file_dir read_file1_pre read_file2_pre target_isize isize_sd repeat_report library_name sample_name read_group_name frag_type(MP|PE|MP_NOPE)" 
	exit(1)

working_dir = sys.argv[1]
reference = sys.argv[2]
read_file_dir = sys.argv[3]
read_file1 = sys.argv[4]
read_file2 = sys.argv[5]
target_isize = sys.argv[6]
isize_sd = sys.argv[7]
repeat_report = sys.argv[8]
library_name = sys.argv[9]
sample_name = sys.argv[10]
read_group_name = sys.argv[11]
frag_type = sys.argv[12]

if (not frag_type in ['MP', 'PE', 'MP_NOPE']):
	print "valid frag types are MP, PE, and MP_NOPE"
	exit(1)

submit_file_dir = "/l2/users/whelanch/gene_rearrange/submit_files"

def prep_submit_file(submit_file_dir, submit_file_name, read_group):
       template = open(submit_file_dir + "/" + submit_file_name, "r")
       outfile = open(submit_file_name, "w")
       for line in template.readlines():
               line = re.sub("(Log\s+=\s+)(\S+)", r'\1/tmp/' + read_group + r'\2', line)
               outfile.write(line)               

read1_files = glob.glob(read_file_dir + read_file1 + "*")
read2_files = glob.glob(read_file_dir + read_file2 + "*")

assert(len(read1_files) > 0), "didn't find any files matching " + read_file_dir + read_file1 + "*"
assert(len(read1_files) == len(read2_files)), "found %d read1 files and %d read2 files; quitting!" % (len(read1_files), len(read2_files))
num_chunks = len(read1_files)

jobs = 0

os.chdir(working_dir)
dagfile = open("{0}/novoalign_first_tier.dag".format(working_dir), 'w')

first_alignment_chunk_job = jobs
for i in xrange(1,num_chunks+1):
	os.mkdir("wd{0}".format(i))
	file1_name = "{0}_{1:03d}.fastq.gz".format(read_file1, i)
	file2_name = "{0}_{1:03d}.fastq.gz".format(read_file2, i)
	#shutil.copy(read_file_dir + file1_name, "wd{0}".format(i))
	os.symlink(read_file_dir + file1_name, "wd{0}/{1}".format(i, file1_name))
	#shutil.copy(read_file_dir + file2_name, "wd{0}".format(i)
	os.symlink(read_file_dir + file2_name, "wd{0}/{1}".format(i, file2_name))
	chunk_job = i + jobs - 1
	if frag_type == "MP":
		submit_file = "novoalign_tier1.desc"
	if frag_type == "PE":
		submit_file = "novoalign_tier1_pe.desc"
	if frag_type == "MP_NOPE":
		submit_file = "novoalign_tier1_mp_no_pe.desc"
	prep_submit_file(submit_file_dir, submit_file, read_group_name)	
	dagfile.write("JOB {0} {3}/{2} DIR {3}/wd{1}/\n".format(chunk_job, i, submit_file, working_dir))
	dagfile.write("VARS {0} read_file1=\"{2}/wd{1}/{3}\"\n".format(chunk_job, i, working_dir, file1_name))
        dagfile.write("VARS {0} read_file2=\"{2}/wd{1}/{3}\"\n".format(chunk_job, i, working_dir, file2_name))
	dagfile.write("VARS {0} reference=\"{1}\"\n".format(chunk_job, reference))
	dagfile.write("VARS {0} target_isize=\"{1}\"\n".format(chunk_job, target_isize))
	dagfile.write("VARS {0} isize_sd=\"{1}\"\n".format(chunk_job, isize_sd))
	dagfile.write("VARS {0} repeat_report=\"{1}\"\n".format(chunk_job, repeat_report))
	dagfile.write("VARS {0} library_name=\"{1}\"\n".format(chunk_job, library_name))
	dagfile.write("VARS {0} sample_name=\"{1}\"\n".format(chunk_job, sample_name))
	dagfile.write("VARS {0} read_group_name=\"{1}\"\n".format(chunk_job, read_group_name))
	dagfile.write("RETRY {0} 2\n".format(chunk_job))
	param_file = open("wd{0}/novo_params".format(i), 'w')
	param_file.write("FORMAT=ILM1.8\n")
	param_file.close()

jobs = jobs + num_chunks 

merge_job = jobs
submit_file = "merge_bams.desc"
prep_submit_file(submit_file_dir, submit_file, read_group_name)
dagfile.write("JOB {0} {1}/{2} DIR {1}\n".format(merge_job, working_dir, submit_file))

dagfile.write("VARS {0} output=\"novoalign_tier1\"\n".format(merge_job))
dagfile.write("VARS {0} working_dir=\"{1}\"\n".format(merge_job, working_dir))
dagfile.write("VARS {0} in_file_names=\"novoalign_sorted.bam\"\n".format(merge_job))

jobs = jobs + 1

merge_calfiles_job = jobs
submit_file = "merge_calfiles.desc"
prep_submit_file(submit_file_dir, submit_file, read_group_name)
dagfile.write("JOB {0} {1}/{2} DIR {1}\n".format(merge_calfiles_job, working_dir, submit_file))

dagfile.write("VARS {0} working_dir=\"{1}\"\n".format(merge_calfiles_job, working_dir))

jobs = jobs + 1

clean_sam_job = jobs
submit_file = "clean_sam.desc"
prep_submit_file(submit_file_dir, submit_file, read_group_name)
dagfile.write("JOB {0} {1}/{2} DIR {1}\n".format(clean_sam_job, working_dir, submit_file))

dagfile.write("VARS {0} input=\"{1}/{2}\"\n".format(clean_sam_job, working_dir, "novoalign_tier1_sort.bam"))

jobs = jobs + 1

mark_dups_job = jobs
submit_file = "mark_dups.desc"
prep_submit_file(submit_file_dir, submit_file, read_group_name)
dagfile.write("JOB {0} {1}/{2} DIR {1}\n".format(mark_dups_job, working_dir, submit_file))

dagfile.write("VARS {0} input=\"{1}/{2}\"\n".format(mark_dups_job, working_dir, "novoalign_tier1_sort_clean.bam"))

jobs = jobs + 1

calculate_insert_sizes_job = jobs
if frag_type == "PE":
	submit_file = "calculate_insert_sizes_pe.desc"
else:
	submit_file = "calculate_insert_sizes_mp.desc"
prep_submit_file(submit_file_dir, submit_file, read_group_name)

dagfile.write("JOB {0} {2}/{1} DIR {2}\n".format(calculate_insert_sizes_job, submit_file, working_dir))
dagfile.write("VARS {0} input=\"{1}/{2}\"\n".format(calculate_insert_sizes_job, working_dir, "novoalign_tier1_sort_clean_rmdup.bam"))

jobs = jobs + 1

cleanup_job = jobs
submit_file = "cleanup_working_dirs.desc"
prep_submit_file(submit_file_dir, submit_file, read_group_name)
dagfile.write("JOB {0} {1}/{2} DIR {1}\n".format(cleanup_job, working_dir, submit_file))
dagfile.write("VARS {0} working_dir=\"{1}\"\n".format(cleanup_job, working_dir))

jobs = jobs + 1

extract_discordants_job = jobs
submit_file = "extract_discordant_reads.desc"
prep_submit_file(submit_file_dir, submit_file, read_group_name)
dagfile.write("JOB {0} {1}/{2} DIR {1}\n".format(extract_discordants_job, working_dir, submit_file))
dagfile.write("VARS {0} bamfile=\"{1}/{2}\"\n".format(extract_discordants_job, working_dir, "novoalign_tier1_sort_clean_mdup_nsort.bam"))

dagfile.write("PARENT ")
for i in xrange(first_alignment_chunk_job,num_chunks+first_alignment_chunk_job):
	dagfile.write("{0} ".format(i))
dagfile.write("CHILD {0} {1}\n".format(merge_job, merge_calfiles_job))
dagfile.write("PARENT {0} CHILD {1}\n".format(merge_job, clean_sam_job))
dagfile.write("PARENT {0} CHILD {1}\n".format(clean_sam_job, mark_dups_job))
dagfile.write("PARENT {0} CHILD {1}\n".format(mark_dups_job, calculate_insert_sizes_job))
dagfile.write("PARENT {0} CHILD {1}\n".format(mark_dups_job, extract_discordants_job))
dagfile.write("PARENT {0} CHILD {1}\n".format(mark_dups_job, cleanup_job))
dagfile.write("PARENT {0} CHILD {1}\n".format(merge_calfiles_job, cleanup_job))

