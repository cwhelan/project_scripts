#!/usr/bin/env python

import os
import sys
import subprocess
import argparse
from multiprocessing import Pool

# import math, datetime
# import rpy2.robjects.lib.ggplot2 as ggplot2
# from rpy2.robjects import *
# from rpy2.robjects.packages import importr
# base = importr('base')
# grdevices = importr('grDevices')

parser = argparse.ArgumentParser()
parser.add_argument("bp_file", help="BED file of breakpoint regions")
parser.add_argument("bp_name", ="Name of this set of breakpoint regions")
parser.add_argument("feature_file", help="BED file of features to test against")
parser.add_argument("feature_name", help="name of the features being tested against")
parser.add_argument("genome", help="FAI file of the genome")
parser.add_argument("--gaps", help="BED file of gaps in the genome")
parser.add_argument("--process_iteration", help="one of bp_hits|feature_hits")
parser.add_argument("--iteration_number", help="iteration number")

args = parser.parse_args()

bp_file = parser.bp_file
bp_name = parser.bp_name
feature_file = parser.feature_file
feature_name = parser.feature_name
genome = parser.genome
gaps = parse.gaps
process_iteration = parser.process_iteration
iteration_number = parser.iteration_number

script_path = os.path.dirname( os.path.realpath( __file__ ) )

def process_iteration(i, fun, gaps, bp_file, genome, feature_file):
    if gaps != None:
        shufflep = subprocess.Popen(['shuffleBed', '-excl', gaps, '-i', bp_file, '-g', genome, '-chrom'], stdout=subprocess.PIPE)
    else:
        shufflep = subprocess.Popen(['shuffleBed', '-i', bp_file, '-g', genome, '-chrom'], stdout=subprocess.PIPE)
    sortp = subprocess.Popen('sort -k1,1 -k2,2n', stdin=shufflep.stdout, stdout=subprocess.PIPE)
    shuffled = sortp.communicate()[0]
    hits = fun(feature_file, 'stdin', input=shuffled)
    #    out_file.write(str(i) + "\t" + str(hits) + "\n")
    return (i, hits)

# doesn't support passing gaps yet
def wrap_process_iteration(arg_tuple):
    return subprocess.Popen(['condor_run', 'python', script_path + 'overlap_bps_with_feature.py', arg_tuple[3], "none", arg_tuple[5], "none", arg_tuple[4], "--process_iteration", arg_tuple[1], "--iteration_number", str(arg_tuple[0])]).communicate[0]

def compute_bases_overlapped(feature_file, bp_file, input=None):
    if bp_file == 'stdin':
        random_bp_hits = subprocess.Popen(['intersectBed', '-a', feature_file, '-b', bp_file, '-wo',], stdin=subprocess.PIPE, stdout=subprocess.PIPE).communicate(input)[0]
    else:
        random_bp_hits = subprocess.Popen(['intersectBed', '-a', feature_file, '-b', bp_file, '-wo',], stdout=subprocess.PIPE).communicate()[0]
    bases_overlapped = 0
    for line in random_bp_hits.split("\n"):        
        fields = line.split("\t")
        if len(fields) > 6:
            bases_overlapped += int(fields[len(fields) - 1])
    return bases_overlapped

def compute_bp_hits(feature_file, bp_file, input=None):
    if bp_file == 'stdin':
        bp_hits = subprocess.Popen(['intersectBed', '-b', feature_file, '-a', bp_file, '-wa', '-u', '-sorted'], stdin=subprocess.PIPE, stdout=subprocess.PIPE).communicate(input)[0]
    else:
        bp_hits = subprocess.Popen(['intersectBed', '-b', feature_file, '-a', bp_file, '-wa', '-u', '-sorted'], stdout=subprocess.PIPE).communicate()[0]
    hits = 0
    current_bp = ""
    bps = set()
    for line in bp_hits.split("\n"):
        fields = line.split("\t")
        if len(fields) >= 6:
            bp = "\t".join(fields[3:6])
            bps.add(bp)
    return len(bps)

def compute_feature_hits(feature_file, bp_file, input=None):
    if bp_file == 'stdin':
        bp_hits = subprocess.Popen(['intersectBed', '-a', feature_file, '-b', bp_file, '-wa', '-u', '-sorted'], stdin=subprocess.PIPE, stdout=subprocess.PIPE).communicate(input)[0]
    else:
        bp_hits = subprocess.Popen(['intersectBed', '-a', feature_file, '-b', bp_file, '-wa', '-u', '-sorted'], stdout=subprocess.PIPE).communicate()[0]
    hits = len(bp_hits.split("\n")) - 1
    return hits

iterations = 10000
p = Pool(4)

if process_iteration != None:

    # compute the real number of bps that overlap a feature
    num_real_bp_hits = compute_bp_hits(feature_file, bp_file)
    print "num bps that overlap features: {0}".format(num_real_bp_hits)

    # compute the real number of features that overlap a bp region
    num_real_feature_hits = compute_feature_hits(feature_file, bp_file)
    print "num features that overlap bps: {0}".format(num_real_feature_hits)

    # compute the real number of bases overlapped
    real_bases_overlapped = compute_bases_overlapped(feature_file, bp_file)
    print "num bases in bps that overlap features: {0}".format(real_bases_overlapped)

    # use shuffleBed to permute bp locations
    print "Computing permutations for BP hits"
    bps_with_hits_file = open('bps_with_hits.txt', 'w')
    num_random_bp_hits = p.map(wrap_process_iteration, [(i, "bp_hits", gaps, bp_file, genome, feature_file) for i in xrange(iterations)])
    for ip in num_random_bp_hits:
        bps_with_hits_file.write(str(ip[0]) + "\t" + str(ip[1]) + "\n")
    bps_with_hits_file.close()                             

    print "Computing permutations for feature hits"
    feature_hits_file = open('feature_hits.txt', 'w')
    num_random_feature_hits = p.map(wrap_process_iteration, [(i, "feature_hits", gaps, bp_file, genome, feature_file) for i in xrange(iterations)])
    for ip in num_random_feature_hits:
        feature_hits_file.write(str(ip[0]) + "\t" + str(ip[1]) + "\n")
    feature_hits_file.close()                            

    print "Computing permutations for bases overlapped"
    bases_overlap_file = open('bases_overlapped.txt', 'w')
    num_bases_overlapped = p.map(wrap_process_iteration, [(i, compute_bases_overlapped, gaps, bp_file, genome, feature_file) for i in xrange(iterations)])
    for ip in num_bases_overlapped:
        bases_overlap_file.write(str(ip[0]) + "\t" + str(ip[1]) + "\n")
    bases_overlap_file.close()                            

    print "Computing shifts"
    # shift the regions over and compute overlaps
    shifted_bp_hits_file = open('shifted_bp_hits.txt', 'w')
    shifted_feature_hits_file = open('shifted_feature_hits.txt', 'w')
    shifted_bases_overlapped_file = open('shifted_bases_overlapped.txt', 'w')
    for i in xrange(-1000000,1000000,10000):
        shifted_regions = subprocess.Popen(['slopBed', '-i', bp_file, '-g', genome, '-l', str(i), '-r', str(-1 * i)], stdout=subprocess.PIPE).communicate()[0]
        bp_hits = compute_bp_hits(feature_file, 'stdin', shifted_regions)
        shifted_bp_hits_file.write(str(i) + "\t" + str(bp_hits) + "\n")
        feature_hits = compute_feature_hits(feature_file, 'stdin', shifted_regions)
        shifted_feature_hits_file.write(str(i) + "\t" + str(feature_hits) + "\n")
        bases_overlapped = compute_bases_overlapped(feature_file, 'stdin', shifted_regions)
        shifted_bases_overlapped_file.write(str(i) + "\t" + str(bases_overlapped) + "\n")
    shifted_bp_hits_file.close()
    shifted_feature_hits_file.close()
    shifted_bases_overlapped_file.close()

    # plot the results
    subprocess.call(map(str, ['/g/whelanch/software/bin/Rscript', script_path + '/plot_results.R', feature_name, bp_name, num_real_bp_hits, num_real_feature_hits, real_bases_overlapped]))

elif process_iteration == "feature_hits":
    print process_iteration(iteration_number, compute_feature_hits, gaps, bp_file, genome)
elif process_iteration == "bp_hits":
    print process_iteration(iteration_number, compute_bp_hits, gaps, bp_file, genome)
else:
    sys.stderr.write("Bad process_iteration commmand!")
