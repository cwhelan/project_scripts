#!/usr/bin/env python

import sys
import os

species_vals = {}
feature_names = set()

for infile in sys.argv[1:len(sys.argv)]:
    infile_basename = os.path.basename(infile)
    species = infile_basename.split("_")[0]
    features_for_species = None
    me_rate_for_species = None
    for line in open(infile, "r"):
        if line.startswith("General"):
            features_for_species = line.rstrip().split("\t")[1:len(line.rstrip().split("\t"))]
            for feature_name in features_for_species:
                feature_names.add(feature_name)
        if line.startswith("Me_rate"):
            me_rate_for_species = line.rstrip().split("\t")[1:len(line.rstrip().split("\t"))]
    assert(len(features_for_species) == len(me_rate_for_species))
    species_map = {}
    for i in xrange(0,len(features_for_species)):
        species_map[features_for_species[i]] = me_rate_for_species[i]
    species_vals[species] = species_map

print "\t".join([""] + species_vals.keys())
for feature_name in feature_names:
    sys.stdout.write(feature_name + "\t")
    for species in species_vals.keys():
        species_map = species_vals[species]
        sys.stdout.write("\t")
        if feature_name in species_map.keys():
            sys.stdout.write(species_map[feature_name])
        else:
            sys.stdout.write("NA")
    sys.stdout.write("\n")
    


    
