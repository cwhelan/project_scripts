#!/usr/bin/env python

import sys
import subprocess

bp_regions_by_file = {"cluster_regions_sorted.bed3.bed": "All Clusters",
                      "cluster_regions_functional_sorted.bed3.bed": "Functional Clusters",
                      "cluster_regions_structural_sorted.bed3.bed": "Structural Clusters",
                      "all_breakpoints_sorted.bed3.bed": "All Breakpoints",
                      }

features_by_file = {
    "ARS-S.cerevisiae.bed.gz.sorted.bed3.bed.gz": "ARS-S.cerevisiae",
    "Alu.bed.gz.sorted.bed3.bed.gz": "Alu",
    "Chi-like.bed.gz.sorted.bed3.bed.gz": "Chi-like",
    "EukaryotesReplicationOrigin.bed.gz.sorted.bed3.bed.gz": "Eukaryotes Replication Origin",
    "HumanMinisatelliteCore.bed.gz.sorted.bed3.bed.gz": "Human Minisatellite Core",
    "LINE.bed.gz.sorted.bed3.bed.gz" : "LINE",
    "PurI.bed.gz.sorted.bed3.bed.gz" : "PurI",
    "PutativeTripleHelices.bed.gz.sorted.bed3.bed.gz" : "Putative Triple Helices",
    "PyrimidineTrait.bed.gz.sorted.bed3.bed.gz" : "Pyrimidine Trait",
    "SINE.bed.gz.sorted.bed3.bed.gz" : "SINE",
    "SatelliteIIIcore.bed.gz.sorted.bed3.bed.gz" : "SatelliteIII core",
    "TopoII_1.bed.gz.sorted.bed3.bed.gz" : "TopoII_1",
    "TopoII_2.bed.gz.sorted.bed3.bed.gz" : "TopoII_2",
    "TranslinMajor.bed.gz.sorted.bed3.bed.gz" : "TranslinMajor",
    "TranslinMajorN.bed.gz.sorted.bed3.bed.gz" : "TranslinMajorN",
    "TranslinMajorNN.bed.gz.sorted.bed3.bed.gz" : "TranslinMajorNN",
    "TranslinMinor.bed.gz.sorted.bed3.bed.gz" : "TranslinMinor",
    "TranslinMinorN.bed.gz.sorted.bed3.bed.gz" : "TranslinMinorN",
    "TranslinMinorNN.bed.gz.sorted.bed3.bed.gz" : "TranslinMinorNN",
    "TranslinMinorNNN.bed.gz.sorted.bed3.bed.gz" : "TranslinMinorNNN",
    "TranslinMinorNNNN.bed.gz.sorted.bed3.bed.gz" : "TranslinMinorNNNN",
    "VDJRSS_7.bed.gz.sorted.bed3.bed.gz" : "VDJRSS_7",
    "VDJRSS_9.bed.gz.sorted.bed3.bed.gz" : "VDJRSS_9",
    "segmental_duplications.bed.gz.sorted.bed3.bed.gz" : "Segmental Duplications",
    "wgEncodeRegDnaseClustered.bed.gz.sorted.bed3.bed.gz" : "ENCODE DNAse Regions"}

print "all: \"
for bp_region in bp_regions_by_file.keys():
    for feature in features_by_file.keys():
        bp_name = bp_regions_by_file.get(bp_region)
        feature_name = features_by_file.get(feature)
        combo_name = bp_name.replace(' ', '_') + "_to_" + feature_name.replace(' ', '_')
        print "\t{0}/{0}.pdf \\".format(combo_name, bp_region, feature)
print

for bp_region in bp_regions_by_file.keys():
    for feature in features_by_file.keys():
        bp_name = bp_regions_by_file.get(bp_region)
        feature_name = features_by_file.get(feature)
        combo_name = bp_name.replace(' ', '_') + "_to_" + feature_name.replace(' ', '_')
        print "{0}/{0}.pdf: {1} {2}".format(combo_name, bp_region, feature)
        print "\tpython overlap_bps_with_feature.py {0} '{1}' {2} '{3}' hg18.fa.fai --permutations 100000 --iteration_chunk_size 100 --cores 100".format(bp_region, bp_name, feature, feature_name)
        print 

