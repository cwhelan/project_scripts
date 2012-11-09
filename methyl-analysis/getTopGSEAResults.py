#!/usr/bin/env python

import glob
import itertools
import subprocess

for gs in ['c2.all.v3.1.symbols.gmt', 'c5.all.v3.1.symbols.gmt', 'c6.all.v3.1.symbols.gmt']:
    hyper = open(gs + "_results_hyper.txt", 'w')
    hypo = open(gs + "_results_hypo.txt", 'w')
    hypervals = {}
    hypovals = {}
    types = ['normalGenes', 'intermediateGenes', 'cancerGenes', 'normalIntermediateGeneDiff', 'intermediateCancerGeneDiff', 'normalCancerGeneDiff', 'normalPromoters', 'intermediatePromoters', 'cancerPromoters', 'normalIntermediatePromoters', 'intermediateCancerPromoters', 'normalCancerPromoters']
    for s in types:
        d1 = s + ".rnk.gsea.out"
        d2 = glob.glob(d1 + "/" + s + "_" + gs + "*")[0]
        posfile = glob.glob(d2 + "/" + "gsea_report_for_na_pos_*.xls")[0]
        pos = subprocess.check_output("cat " + posfile + "| cut -f1,6,8 | awk '$3 < 0.25' | sort -k2,2nr | cut -f1", shell=True)
        negfile = glob.glob(d2 + "/" + "gsea_report_for_na_neg_*.xls")[0]
        neg = subprocess.check_output("cat " + negfile + "| cut -f1,6,8 | awk '$2 != 1.0 && $3 < 0.25' | sort -k2,2n | cut -f1", shell=True)
        hypervals[s] = pos.split("\n")
        hypovals[s] = neg.split("\n")
    for x in itertools.izip_longest(*[hypervals[s] for s in types], fillvalue=""):
        hyper.write("\t".join(str(i) for i in x) + "\n")
    for x in itertools.izip_longest(*[hypovals[s] for s in types], fillvalue=""):
        hypo.write("\t".join(str(i) for i in x) + "\n")        

