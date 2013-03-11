#!/usr/bin/env python

import os
import sys
import subprocess
from multiprocessing import Pool

# import math, datetime
# import rpy2.robjects.lib.ggplot2 as ggplot2
# from rpy2.robjects import *
# from rpy2.robjects.packages import importr
# base = importr('base')
# grdevices = importr('grDevices')

bp_file = sys.argv[1]
bp_name=sys.argv[2]
feature_file = sys.argv[3]
feature_name = sys.argv[4]

genome = sys.argv[5]
gaps = None
if len(sys.argv) == 7:
    gaps = sys.argv[6]

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

def wrap_process_iteration(arg_tuple):
    return process_iteration(*arg_tuple)

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
num_random_bp_hits = p.map(wrap_process_iteration, [(i, compute_bp_hits, gaps, bp_file, genome, feature_file) for i in xrange(iterations)])
#for i in xrange(iterations):    
    # process_iteration(i, compute_bp_hits, bps_with_hits_file, gaps, bp_file, genome, feature_file, num_random_bp_hits)
    # shufflep = subprocess.Popen(['shuffleBed', '-excl', gaps, '-i', bp_file, '-g', genome, '-chrom'], stdout=subprocess.PIPE)
    # hits = compute_bp_hits(feature_file, 'stdin', shufflep.communicate()[0])
    # num_random_bp_hits.append(hits)
    # bps_with_hits_file.write(str(i) + "\t" + str(hits) + "\n")
for ip in num_random_bp_hits:
    bps_with_hits_file.write(str(ip[0]) + "\t" + str(ip[1]) + "\n")
bps_with_hits_file.close()                             

print "Computing permutations for feature hits"
feature_hits_file = open('feature_hits.txt', 'w')
num_random_feature_hits = p.map(wrap_process_iteration, [(i, compute_feature_hits, gaps, bp_file, genome, feature_file) for i in xrange(iterations)])
#for i in xrange(iterations):
    # process_iteration(i, compute_feature_hits, feature_hits_file, gaps, bp_file, genome, feature_file, num_random_feature_hits)
    # shufflep = subprocess.Popen(['shuffleBed', '-excl', gaps, '-i', bp_file, '-g', genome, '-chrom'], stdout=subprocess.PIPE)
    # hits = compute_feature_hits(feature_file, 'stdin', shufflep.communicate()[0])
    # num_random_feature_hits.append(hits)
    # feature_hits_file.write(str(i) + "\t" + str(hits) + "\n")
for ip in num_random_feature_hits:
    feature_hits_file.write(str(ip[0]) + "\t" + str(ip[1]) + "\n")
feature_hits_file.close()                            

print "Computing permutations for bases overlapped"
bases_overlap_file = open('bases_overlapped.txt', 'w')
num_bases_overlapped = p.map(wrap_process_iteration, [(i, compute_bases_overlapped, gaps, bp_file, genome, feature_file) for i in xrange(iterations)])
# for i in xrange(iterations):
    # process_iteration(i, compute_bases_overlapped, bases_overlap_file, gaps, bp_file, genome, feature_file, num_bases_overlapped)
    # shufflep = subprocess.Popen(['shuffleBed', '-excl', gaps, '-i', bp_file, '-g', genome, '-chrom'], stdout=subprocess.PIPE)
    # bases_overlapped = compute_bases_overlapped(feature_file, 'stdin', shufflep.communicate()[0])
    # num_bases_overlapped.append(bases_overlapped)
    # bases_overlap_file.write(str(i) + "\t" + str(bases_overlapped) + "\n")
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
