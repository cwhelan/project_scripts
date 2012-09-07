#!/bin/env python

# this script takes the text file of lava insertions generated by alan harris and converts it
# into a flat file importable into R

# the file comes in this format:

# lines that start with a '#' contain information about a repeat location:
#Repeat_ID	Repeat_Scaffold	Repeat_Start	Repeat_Stop

# NB it looks like the data is actually in the format:
#Repeat_Scaffold	Repeat_Start	Repeat_Stop    Repeat_id   "OVERLAP"

# those are followed by one or more lines like this
#Gene_Name	Gene_Type	Transcript_ID	Gene_ID	Protiein_ID	Scaffold	Start	Stop	Strand	Gene_Location	CDS

# this script puts all that info on one line:
#Repeat_ID	Repeat_Scaffold	Repeat_Start	Repeat_Stop Gene_Name	Gene_Type	Transcript_ID	Gene_ID	Protiein_ID	Scaffold	Start	Stop	Strand	Gene_Location	CDS

import sys

input_file = open(sys.argv[1], "r")
output_file = open(sys.argv[2], "w")

output_file.write("Repeat_ID\tRepeat_Scaffold\tRepeat_Start\tRepeat_Stop\tGene_Name\tGene_Type\tTranscript_ID\tGene_ID\tProtiein_ID\tScaffold\tStart\tStop\tStrand\tGene_Location\tCDS\n")

current_lava_line = "NO_CURRENT_LAVA_LINE"
for line in input_file:
    if line.startswith("#"):
        current_lava_fields = line.rstrip()[1:].split("\t")
        current_lava_line = "\t".join([current_lava_fields[3], current_lava_fields[0], current_lava_fields[1], current_lava_fields[2]])
    else:
        output_file.write(current_lava_line + "\t" + line.rstrip() + "\n")
