#!/usr/bin/env python

import argparse
<<<<<<< HEAD
import sys
=======
>>>>>>> 5fabbd20142fbc40b1f140f5bbf6daa7df8eac93

parser = argparse.ArgumentParser()
parser.add_argument("event_name", help="name of the event from the bedpe file")
parser.add_argument("supporting_reads_file", help="filename of the breakdancer supporting reads file")

args = parser.parse_args()

event_name = args.event_name
in_event = False
<<<<<<< HEAD
for line in open(args.supporting_reads_file):
    if line.startswith("track"):
        if in_event:
            sys.exit(0)
=======
for line in sys.open(args.supporting_reads_file):
    if line.startswith("track"):
        if in_event:
            exit
>>>>>>> 5fabbd20142fbc40b1f140f5bbf6daa7df8eac93
        else:
            this_event_name = line.split()[1].split("=")[1]
            if (this_event_name == event_name):
                in_event = True
<<<<<<< HEAD
                sys.stdout.write(line)
	    else:
		in_event = False
    else:
        if in_event:
            sys.stdout.write(line)
=======
                print line
    else:
        if in_event:
            print line
        
>>>>>>> 5fabbd20142fbc40b1f140f5bbf6daa7df8eac93
