#!/usr/local/bin/python                                                                                                                           

import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def make_rc_record(record):
    """Returns a new SeqRecord with the reverse complement sequence."""
    record.letter_annotations['phred_quality'].reverse()
    return SeqRecord(seq = record.seq.reverse_complement(), \
                     id = "rc_" + record.id, \
                     letter_annotations = record.letter_annotations,
		     description = "")

trimmed_reads = (make_rc_record(rec) for rec in \
		       SeqIO.parse(sys.argv[1], "fastq"))
count = SeqIO.write(trimmed_reads, sys.argv[2], "fastq")
print "Saved %i reads" % count


