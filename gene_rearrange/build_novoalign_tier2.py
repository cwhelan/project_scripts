#!/usr/bin/env python

import sys
import os
import shutil
import subprocess
import time
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("working_dir", help="Working directory")
parser.add_argument("reference", help="Novoalign reference file")
parser.add_argument("read_file1", help="fastq file of read 1 of read pairs to be realigned")
parser.add_argument("read_file2", help="fastq file of read 2 of read pairs to be realigned")
parser.add_argument("target_isize", help="library insert size", type=int)
parser.add_argument("isize_sd", help="library insert size standard deviation", type=int)
parser.add_argument("tier1_stat_file", help="statistics file for tier 1 alignment")
parser.add_argument("threshold", help="Novoalign alignment threshold", type=int)
parser.add_argument("frag_type", help="Fragment type (MP|PE|MP_NOPE)")
args = parser.parse_args()

working_dir = args.working_dir
reference = args.reference
read_file1 = args.read_file1
read_file2 = args.read_file2
target_isize = str(args.target_isize)
isize_sd = str(args.isize_sd)
tier1_stat_file = args.tier1_stat_file
threshold = str(args.threshold)
frag_type = args.frag_type

chunk_size = 100000

def split_file(name, prefix, chunk_size):
        split_cmd = "split -a 4 -d -l {2} {0} {1}.".format(name, prefix, chunk_size)
        subprocess.call(split_cmd, shell=True)

wcout = subprocess.Popen("wc -l {0}".format(read_file1), shell=True, stdout=subprocess.PIPE).communicate()[0]
file_length = int(wcout.split()[0])
num_chunks = (int(file_length) / chunk_size) + 1

os.chdir(working_dir)
dagfile = open("{0}/novoalign_second_tier.dag".format(working_dir), 'w')
file1_base_name = os.path.basename(read_file1)
split_file(read_file1, file1_base_name, chunk_size)
file2_base_name = os.path.basename(read_file2)
split_file(read_file2, file2_base_name, chunk_size)

jobs = 0

first_alignment_chunk_job = jobs
for i in xrange(jobs,num_chunks+jobs):
	os.mkdir("wd{0}".format(i))
	shutil.move("{0}.{1:04d}".format(file1_base_name, i), "wd{0}".format(i))
	shutil.move("{0}.{1:04d}".format(file2_base_name, i), "wd{0}".format(i))

	if frag_type == "MP":
		submit_file = "novoalign_tier2.desc"
	if frag_type == "PE":
		submit_file = "novoalign_tier2_pe.desc"
	if frag_type == "MP_NOPE":
		submit_file = "novoalign_tier2_mp_no_pe.desc"

	dagfile.write("JOB {0} /l2/users/whelanch/gene_rearrange/submit_files/{1} DIR wd{0}/\n".format(i, submit_file))

	dagfile.write("VARS {0} read_file1=\"{1}/wd{0}/{2}.{0:04d}\"\n".format(i, working_dir, file1_base_name))
        dagfile.write("VARS {0} read_file2=\"{1}/wd{0}/{2}.{0:04d}\"\n".format(i, working_dir, file2_base_name))
	dagfile.write("VARS {0} reference=\"{1}\"\n".format(i, reference))
	dagfile.write("VARS {0} target_isize=\"{1}\"\n".format(i, target_isize))
	dagfile.write("VARS {0} isize_sd=\"{1}\"\n".format(i, isize_sd))
	dagfile.write("VARS {0} threshold=\"{1}\"\n".format(i, threshold))
	dagfile.write("RETRY {0} 2\n".format(i))

jobs = jobs + num_chunks

merge_job = jobs
dagfile.write("JOB {0} /l2/users/whelanch/gene_rearrange/submit_files/merge_bams.desc DIR {1}\n".format(merge_job, working_dir))
dagfile.write("VARS {0} output=\"novoalign_tier2_nsort\"\n".format(merge_job))
dagfile.write("VARS {0} working_dir=\"{1}\"\n".format(merge_job, working_dir))
dagfile.write("VARS {0} in_file_names=\"novoalign_sorted.bam\"\n".format(merge_job))

jobs = jobs + 1

clean_sam_job = jobs
dagfile.write("JOB {0} /l2/users/whelanch/gene_rearrange/submit_files/clean_sam.desc DIR {1}\n".format(clean_sam_job, working_dir))
dagfile.write("VARS {0} input=\"{1}/{2}\"\n".format(clean_sam_job, working_dir, "novoalign_tier2_nsort.bam"))

jobs = jobs + 1

hydra_submit_file = "hydra_merged.desc"
if frag_type == "PE":
	hydra_submit_file = "hydra_merged_pe.desc"

hydra_job = jobs
dagfile.write("JOB {0} /l2/users/whelanch/gene_rearrange/submit_files/{1} DIR {2}\n".format(hydra_job, hydra_submit_file, working_dir))
dagfile.write("VARS {0} input=\"{1}/{2}\"\n".format(hydra_job, working_dir, "novoalign_tier2_nsort.bam"))
dagfile.write("VARS {0} tier1_stat_file=\"{1}/{2}\"\n".format(hydra_job, working_dir, "novoalign_tier1_sort_clean_rmdup_insert_stats.txt"))

jobs = jobs + 1

#cleanup_job = jobs
#dagfile.write("JOB {0} /l2/users/whelanch/gene_rearrange/submit_files/cleanup_working_dirs.desc DIR {1}\n".format(cleanup_job, working_dir))
#dagfile.write("VARS {0} working_dir=\"{1}\"\n".format(cleanup_job, working_dir))

#jobs = jobs + 1

dagfile.write("PARENT ")
for i in xrange(first_alignment_chunk_job,num_chunks+first_alignment_chunk_job):
	dagfile.write("{0} ".format(i))
dagfile.write("CHILD {0}\n".format(num_chunks+first_alignment_chunk_job))

dagfile.write("PARENT {0} CHILD {1}\n".format(merge_job, clean_sam_job))
dagfile.write("PARENT {0} CHILD {1}\n".format(clean_sam_job, hydra_job))
#dagfile.write("PARENT {0} CHILD {1}\n".format(hydra_job, cleanup_job))
