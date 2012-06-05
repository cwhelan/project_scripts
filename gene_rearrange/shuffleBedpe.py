#!/usr/bin/env python

import sys
import os
import pybedtools

if len(sys.argv) != 4:
    print "Usage: shuffleBedpe.py bedpe_file genome_file exclusion_region_file"
    sys.exit()

bedpe_file = sys.argv[1]
genome_file = sys.argv[2]
excl_file = sys.argv[3]

bedpe_entries = pybedtools.BedTool(bedpe_file)

num_fields = len(bedpe_entries[0].fields)

def outer_span_feature(feature):
    #print feature
    return pybedtools.Interval(feature.chrom, feature.start, long(feature.fields[5]), feature.fields[6], feature.fields[7], '+', [feature.fields[1], feature.fields[2], feature.fields[4], feature.fields[5]] + feature.fields[8:])

def print_bedpe_feature_from_ospan(feature):
    #print feature
    
    ispan1_len = long(feature.fields[7]) - long(feature.fields[6])
    ispan2_len = long(feature.fields[9]) - long(feature.fields[8])
    print "\t".join(str(x) for x in [feature.chrom, feature.start, feature.start + ispan1_len, feature.chrom, feature.end - ispan2_len, feature.end, feature.name, feature.score] + feature.fields[10:])

def feature_pair_from_bedpe(feature_list):
    for feature in feature_list:
        #print feature
        #print feature.fields
        f1 = pybedtools.Interval(feature.chrom, feature.start, feature.end, feature.fields[6] + "-1", feature.fields[7], feature.fields[8], feature.fields[10:])
        f2 = pybedtools.Interval(feature.fields[3], long(feature.fields[4]), long(feature.fields[5]), feature.fields[6] + "-2", feature.fields[7], feature.fields[9], feature.fields[10:])
        yield f1
        yield f2

def print_bedpe_feature_from_end_pairs(f1, f2):
    print "\t".join(str(x) for x in [f1.chrom, f1.start, f1.end, f2.chrom, f2.start, f2.end, f1.name.split("-")[0], f1.score, f1.strand, f2.strand] + f1.fields[6:])

intra_chrom_ospans = pybedtools.BedTool(outer_span_feature(f) for f in bedpe_entries if f.fields[0] == f.fields[3])
shuffled_intra_chrom_ospans = intra_chrom_ospans.shuffle(g=genome_file, chrom=True, excl=excl_file)

inter_chrom_ends = pybedtools.BedTool(feature_pair_from_bedpe((f for f in bedpe_entries if f.fields[0] != f.fields[3])))
shuffled_inter_chrom_ends = inter_chrom_ends.shuffle(g=genome_file, chrom=True, excl=excl_file)

for feature in shuffled_intra_chrom_ospans:
    print_bedpe_feature_from_ospan(feature)

shuffled_inter_chrom_ends_by_name = {}
for feature in shuffled_inter_chrom_ends:
    if feature.name.split("-")[0] in shuffled_inter_chrom_ends_by_name:
        other = shuffled_inter_chrom_ends_by_name[feature.name.split("-")[0]]
        if feature.name.split("-")[1] == "1":
            print_bedpe_feature_from_end_pairs(feature, other)
        else:
            print_bedpe_feature_from_end_pairs(other, feature)
        del shuffled_inter_chrom_ends_by_name[feature.name.split("-")[0]]
    else:
        shuffled_inter_chrom_ends_by_name[feature.name.split("-")[0]] = feature

        


