#!/user/bin/env python

import sys
import fileinput

selected_id_file = sys.argv[1]
bd_file = sys.argv[2]

selected_ids = []
for line in open(selected_id_file).readlines():
    selected_ids.append( line.strip())

#print selected_ids

in_selected_bp = False
for line in fileinput.input(bd_file):
    line = line.strip()
    if line.startswith("track"):
        name = line.split()[1].split("=")[1]
        #print name
        if name in selected_ids:
            print line
            in_selected_bp = True
        else:
            in_selected_bp = False
    else:
        if in_selected_bp:
            print line
