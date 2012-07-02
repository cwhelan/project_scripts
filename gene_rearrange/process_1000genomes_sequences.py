#!/usr/bin/env python

import os
import sys
import subprocess

sequence_index = open(sys.argv[1], 'r')

scripts_dir = '/l2/users/whelanch/project_scripts/gene_rearrange/'
index_dir = '/l2/users/whelanch/genome_refs/10kg/hg19/'
thousand_genomes_data = '/l2/users/whelanch/1000genomes/data/pilot_data/'

current_dir = os.getcwd()

master_dag = open('master.dag', 'w')

for line in sequence_index:
    fields = line.split("\t")
    seq_file = fields[0]
    paired_file = fields[20]
    if not seq_file.find("_1") > -1 and paired_file.find("_2") > -1:
        continue
    read_group = fields[2] # run_id
    sample = fields[8]
    library = fields[14]
    isize = fields[17]
    os.mkdir(read_group)
    os.chdir(read_group)
    subprocess.Popen([scripts_dir + 'build_novoalign_tier1.py',
                      current_dir,
                      thousand_genomes_data + seq_file,
                      thousand_genomes_data + paired_file,
                      isize,
                      str(float(isize) * .15),
                      'Random',
                      library,
                      sample,
                      read_group,
                      'PE',
                      'STDFQ']).communicate()[0]
    subprocess.Popen(['condor_submit_dag', '-no_submit', 'novoalign_first_tier.dag'])
    master_dag.write("JOB\t" + read_group + "\t" + current_dir + "/" + read_group + "/novoalign_first_tier.dag.condor.sub")
    os.chdir(current_dir)

master_dag.close()
    
    
    
