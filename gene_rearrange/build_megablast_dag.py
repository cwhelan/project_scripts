#!/usr/local/bin/python

import sys
import os
import shutil
import subprocess

os.chdir(sys.argv[1])
dagfile = open("{0}/megablast.dag".format(sys.argv[1]), 'w')
split_cmd = "split -a 3 -d -l 2000 {0} {0}.".format(sys.argv[2])
subprocess.Popen(split_cmd, shell=True)
for i in xrange(0,100):
	os.mkdir("{0}".format(i))
	shutil.move("tier1_mapped_sample.fa.{0:03d}".format(i), "{0}".format(i))
	dagfile.write("JOB {0} /l2/users/whelanch/gene_rearrange/scripts/megablast.desc DIR {1}/{0}/\n".format(i, sys.argv[1]))
