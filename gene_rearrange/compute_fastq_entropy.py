import sys
from math import *

def increment(map, key):
    if key in map:
        map[key] += 1
    else:
        map[key] = 1

def entropy(counts):
    normalized_counts = map(lambda x: float(x)/sum(counts.values()), counts.values())
    entropy = -1 * sum(map(lambda x: x * log(x), normalized_counts))
    return entropy

while True:
    header = sys.stdin.readline()
    if not header:
        break
    seq = sys.stdin.readline()
    sep = sys.stdin.readline()
    qual = sys.stdin.readline()

    uni_counts = {}
    bi_counts = {}
    tri_counts = {}

    for i in range(0, len(seq) - 1):
        increment(uni_counts, seq[i])
        if i > 0:
            increment(bi_counts, seq[i-1:i+1])
        if i > 1:
            increment(tri_counts, seq[i-2:i+1])
    print "\t".join(map(str, [entropy(uni_counts), entropy(bi_counts), entropy(tri_counts)]))
    
        
