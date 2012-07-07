#!/usr/bin/env python

import sys
import os
import shutil
import subprocess
import time
from cStringIO import StringIO

if not len(sys.argv) == 13:
	print "Usage: build_novoalign_tier1.py working_dir reference read_file1 read_file2 target_isize isize_sd repeat_report library_name sample_name read_group_name frag_type(MP|PE|MP_NOPE) qual_format" 
	exit(1)

working_dir = sys.argv[1]
reference = sys.argv[2]
read_file1 = sys.argv[3]
read_file2 = sys.argv[4]
target_isize = sys.argv[5]
isize_sd = sys.argv[6]
repeat_report = sys.argv[7]
library_name = sys.argv[8]
sample_name = sys.argv[9]
read_group_name = sys.argv[10]
frag_type = sys.argv[11]
qual_format = sys.argv[12]

if (not frag_type in ['MP', 'PE', 'MP_NOPE']):
	print "valid frag types are MP, PE, and MP_NOPE"
	exit(1)

chunk_size = 1000000

def split_file(name, prefix, chunk_size):
	if name.endswith("gz"):
		basename = prefix.split(".gz")[0]
		lines_read = 0
		chunk = 0
		p = subprocess.Popen(["zcat",name], 
				     stdout = subprocess.PIPE)
		read_file = StringIO(p.communicate()[0])
		assert p.returncode == 0 
		chunk_filename = "{0}.{1:04}.gz".format(basename,chunk)
		print "opening {0} for writing".format(chunk_filename)
		p = subprocess.Popen("gzip -c > " + chunk_filename, shell=True, stdin=subprocess.PIPE)
		chunk_file = p.stdin
		for line in read_file:
			chunk_file.write(line)
			lines_read = lines_read + 1
			if lines_read % 100000 == 0:
				print "wrote {0} lines".format(lines_read)
			if lines_read % chunk_size == 0:
				chunk_file.close()
				chunk = chunk + 1
				chunk_filename = "{0}.{1:04}.gz".format(basename,chunk)
				print "opening {0} for writing".format(chunk_filename)
				p.communicate()
				p = subprocess.Popen("gzip -c > " + chunk_filename, shell=True, stdin=subprocess.PIPE)
				chunk_file = p.stdin
		p.communicate()
		return lines_read
	else:
		split_cmd = "split -a 4 -d -l {2} {0} {1}.".format(name, prefix, chunk_size)
		print "executing: " + split_cmd
		subprocess.call(split_cmd, shell=True)
		return 1

def prep_submit_file(submit_file_dir, submit_file_name, read_group):
	template = open(submit_file_dir + "/" submit_file_name, "r")
	outfile = open(submit_file_name, "w")
	for line in template.readlines():
		line = re.sub("(Log\s+=\s+)(\S+)", r'\1/tmp/' + read_group + r'/\2', line)
	outfile.write(line)
		

if not read_file1.endswith("gz"):
	wcout = subprocess.Popen("wc -l {0}".format(read_file1), shell=True, stdout=subprocess.PIPE).communicate()[0]
	file_length = int(wcout.split()[0])
	print "Read file length: {0}".format(file_length)

	num_chunks = (int(file_length) / chunk_size) + (1 if file_length % chunk_size > 0 else 0)
	print "Number of chunks: {0}".format(num_chunks)

os.chdir(working_dir)
dagfile = open("{0}/novoalign_first_tier.dag".format(working_dir), 'w')
file1_base_name = os.path.basename(read_file1)
lines_read = split_file(read_file1, file1_base_name, chunk_size)
file2_base_name = os.path.basename(read_file2)
split_file(read_file2, file2_base_name, chunk_size)

gzipped = False
if read_file1.endswith("gz"):
	num_chunks = (int(lines_read) / chunk_size) + (1 if lines_read % chunk_size > 0 else 0)
	gzipped = True

# sleeping a bit to let the file system catch up
#time.sleep(120)

jobs = 0

submit_file_dir = "/l2/users/whelanch/gene_rearrange/submit_files"
prep_submit_file(submit_file_dir, "fastqc.desc", read_group_name)
dagfile.write("JOB {0} {1}{2}  DIR {1}/\n".format(jobs, working_dir, "fastqc.desc"))
dagfile.write("VARS {0} read_file=\"{1}\"\n".format(jobs, read_file1))

jobs = jobs + 1

dagfile.write("JOB {0} {1}{2} DIR {1}/\n".format(jobs, working_dir, "fastqc.desc"))
dagfile.write("VARS {0} read_file=\"{1}\"\n".format(jobs, read_file2))

jobs = jobs + 1
first_alignment_chunk_job = jobs

file1_prefix = file1_base_name
file2_prefix = file2_base_name
if gzipped:
	file1_prefix = file1_base_name.split(".gz")[0]
	file2_prefix = file2_base_name.split(".gz")[0]

for i in xrange(0,num_chunks):
	os.mkdir("wd{0}".format(i))
	if gzipped:
		shutil.move("{0}.{1:04d}.gz".format(file1_prefix, i), "wd{0}".format(i))
		shutil.move("{0}.{1:04d}.gz".format(file2_prefix, i), "wd{0}".format(i))	
	else:
		shutil.move("{0}.{1:04d}".format(file1_prefix, i), "wd{0}".format(i))
		shutil.move("{0}.{1:04d}".format(file2_prefix, i), "wd{0}".format(i))
	chunk_job = i + jobs
	if frag_type == "MP":
		submit_file = "novoalign_tier1.desc"
	if frag_type == "PE":
		submit_file = "novoalign_tier1_pe.desc"
	if frag_type == "MP_NOPE":
		submit_file = "novoalign_tier1_mp_no_pe.desc"
	prep_submit_file(submit_file_dir, submit_file, read_group_name)
	
	dagfile.write("JOB {0} {3}{2} DIR {3}wd{1}/\n".format(chunk_job, i, submit_file, working_dir))
	if gzipped:
		dagfile.write("VARS {0} read_file1=\"{2}/wd{1}/{3}.{1:04d}.gz\"\n".format(chunk_job, i, working_dir, file1_prefix))
		dagfile.write("VARS {0} read_file2=\"{2}/wd{1}/{3}.{1:04d}.gz\"\n".format(chunk_job, i, working_dir, file2_prefix))
	else:
		dagfile.write("VARS {0} read_file1=\"{2}/wd{1}/{3}.{1:04d}\"\n".format(chunk_job, i, working_dir, file1_prefix))
		dagfile.write("VARS {0} read_file2=\"{2}/wd{1}/{3}.{1:04d}\"\n".format(chunk_job, i, working_dir, file2_prefix))
	dagfile.write("VARS {0} reference=\"{1}\"\n".format(chunk_job, reference))
	dagfile.write("VARS {0} target_isize=\"{1}\"\n".format(chunk_job, target_isize))
	dagfile.write("VARS {0} isize_sd=\"{1}\"\n".format(chunk_job, isize_sd))
	dagfile.write("VARS {0} repeat_report=\"{1}\"\n".format(chunk_job, repeat_report))
	dagfile.write("VARS {0} library_name=\"{1}\"\n".format(chunk_job, library_name))
	dagfile.write("VARS {0} sample_name=\"{1}\"\n".format(chunk_job, sample_name))
	dagfile.write("VARS {0} read_group_name=\"{1}\"\n".format(chunk_job, read_group_name))
	
	dagfile.write("RETRY {0} 2\n".format(chunk_job))
	param_file = open("wd{0}/novo_params".format(i), 'w')
	param_file.write("FORMAT={0}\n".format(qual_format))
	param_file.close()


jobs = jobs + num_chunks

merge_job = jobs
submit_file = "merge_bams.desc"
prep_submit_file(submit_file_dir, submit_file, read_group_name)
dagfile.write("JOB {0} {1}{2} DIR {1}\n".format(merge_job, working_dir, submit_file))
dagfile.write("VARS {0} output=\"novoalign_tier1\"\n".format(merge_job))
dagfile.write("VARS {0} working_dir=\"{1}\"\n".format(merge_job, working_dir))
dagfile.write("VARS {0} in_file_names=\"novoalign_sorted.bam\"\n".format(merge_job))

jobs = jobs + 1

merge_calfiles_job = jobs
submit_file = "merge_calfiles.desc"
prep_submit_file(submit_file_dir, submit_file, read_group_name)
dagfile.write("JOB {0} {1}{2} DIR {1}\n".format(merge_calfiles_job, working_dir, submit_file))
dagfile.write("VARS {0} working_dir=\"{1}\"\n".format(merge_calfiles_job, working_dir))

jobs = jobs + 1

clean_sam_job = jobs
submit_file = "clean_sam.desc"
prep_submit_file(submit_file_dir, submit_file, read_group_name)
dagfile.write("JOB {0} {1}{2} DIR {1}\n".format(clean_sam_job, working_dir, submit_file))
dagfile.write("VARS {0} input=\"{1}/{2}\"\n".format(clean_sam_job, working_dir, "novoalign_tier1_sort.bam"))

jobs = jobs + 1

mark_dups_job = jobs
submit_file = "mark_dups.desc"
prep_submit_file(submit_file_dir, submit_file, read_group_name)
dagfile.write("JOB {0} {1}{2} DIR {1}\n".format(mark_dups_job, working_dir, submit_file))
dagfile.write("VARS {0} input=\"{1}/{2}\"\n".format(mark_dups_job, working_dir, "novoalign_tier1_sort_clean.bam"))

jobs = jobs + 1

calculate_insert_sizes_job = jobs
if frag_type == "PE":
	submit_file = "calculate_insert_sizes_pe.desc"
else:
	submit_file = "calculate_insert_sizes_mp.desc"
prep_submit_file(submit_file_dir, submit_file, read_group_name)

dagfile.write("JOB {0} {2}{1} DIR {2}\n".format(calculate_insert_sizes_job, submit_file, working_dir))
dagfile.write("VARS {0} input=\"{1}/{2}\"\n".format(calculate_insert_sizes_job, working_dir, "novoalign_tier1_sort_clean_rmdup.bam"))

jobs = jobs + 1

cleanup_job = jobs
submit_file = "cleanup_working_dirs.desc"
prep_submit_file(submit_file_dir, submit_file, read_group_name)
dagfile.write("JOB {0} {1}{2} DIR {1}\n".format(cleanup_job, working_dir, submit_file))
dagfile.write("VARS {0} working_dir=\"{1}\"\n".format(cleanup_job, working_dir))

jobs = jobs + 1

extract_discordants_job = jobs
submit_file = "extract_discordant_reads.desc"
prep_submit_file(submit_file_dir, submit_file, read_group_name)
dagfile.write("JOB {0} {1}{2} DIR {1}\n".format(extract_discordants_job, working_dir, submit_file))
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


