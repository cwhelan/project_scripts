#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("event_name", help="name of the event from the bedpe file")
parser.add_argument("supporting_reads_file", help="filename of the breakdancer supporting reads file")

args = parser.parse_args()

event_name = args.event_name
in_event = False
for line in sys.open(args.supporting_reads_file):
    if line.startswith("track"):
        if in_event:
            exit
        else:
            this_event_name = line.split()[1].split("=")[1]
            if (this_event_name == event_name):
                in_event = True
                print line
    else:
        if in_event:
            print line
        
