#!/usr/bin/env python

import sys
import os
import pybedtools
import subprocess

if len(sys.argv) != 9:
    print "Usage: annotate_hydra_out.py hydra_file te_file common_deletions_file insert_size output_dir seg_dups cent_tel sample_name"
    sys.exit()


hydra_file = sys.argv[1]
print 'analyzing hydra file: ' + hydra_file
te_file = sys.argv[2]
common_deletions_file = sys.argv[3]
insert_size = int(sys.argv[4])
output_dir = sys.argv[5]
seg_dups_file = sys.argv[6]
cent_tel_file = sys.argv[7]
sample_name = sys.argv[8]

log = open(output_dir + "/annotate.log", "w")
log.write("input file: {0}\n".format(hydra_file))
log.write("sample name: ".format(sample_name))
log.write("te file: {0}\n".format(te_file))
log.write("common deletions file: {0}\n".format(common_deletions_file))
log.write("insert size: {0}\n".format(insert_size))
log.write("segmental duplications file: {0}\n".format(seg_dups_file))
log.write("centromeres and telomeres file: {0}\n".format(cent_tel_file))

#pybedtools.set_tempdir('/l2/users/whelanch/scratch')

common_deletion_overlap_pct = 0.5

def bedpe_lt_length_filter(feature, length):
    if int(feature[4]) - int(feature[2]) < length:
        return True
    return False

def bedpe_gt_length_filter(feature, length):
    if int(feature[4]) - int(feature[2]) > length:
        return True
    return False

def inter_chr_filter(feature):
    if feature[0] == feature[3]:
        return False
    return True

def intra_chr_filter(feature):
    if feature[0] == feature[3]:
        return True
    return False

def score_gte_filter(feature, score):
    if int(feature[7]) >= score:
        return True
    return False

def score_lt_filter(feature, score):
    if int(feature[7]) < score:
        return True
    return False

def expected_orientation_filter(feature, matches):
    expected = False
#    print feature
    if (feature[8] == '+' and feature[9] == '-') or (feature[8] == '-' and feature[9] == '+'):
#        print "found a good feature"
        expected = True

    if matches:
        return expected
    else:
        return not expected

def bedpe_reciprocal_overlap_ends_filter(feature, overlap_pct):
    #print feature
    te_chr = feature[22]
    te_start = int(feature[23])
    te_end = int(feature[24])
    te_length = te_end - te_start
#    match = False
    return overlaps_by(te_chr, te_start, te_end, feature[0], int(feature[1]), int(feature[2]), overlap_pct) or overlaps_by(te_chr, te_start, te_end, feature[3], int(feature[4]), int(feature[5]), overlap_pct)

def bedpe_reciprocal_overlap_ispan_filter(feature, overlap_pct):
 #   print feature
    te_chr = feature[22]
    te_start = int(feature[23])
    te_end = int(feature[24])
    te_length = te_end - te_start
#    match = False
    return overlaps_by(te_chr, te_start, te_end, feature[0], int(feature[2]), int(feature[4]), overlap_pct) 

def overlaps_by(chr1, start1, end1, chr2, start2, end2,  overlap_pct):
#    print [chr1,start1,end1,chr2,start2,end2]
    if chr1 != chr2:
#        print "chrs don't match"
        return False
    max_start = max(start1,start2)
    min_end = min(end1,end2)
#    print "max start %d; min end: %d" % (max_start, min_end)
    overlap_len = min_end - max_start
#    print "overlap_len: %d" % overlap_len
    len1 = end1 - start1
    len2 = end2 - start2
    return float(overlap_len) / len1 >= overlap_pct and float(overlap_len) / len2 >= overlap_pct
    
def write_bed(call, fh):
    for f in call.fields:
        fh.write(str(f) + "\t")
    fh.write("\n")
    
def merge_duplicate_breaks(calls, slop):
    calls_list = []

    for call in calls:
        calls_list.append(call)
    for curr in calls_list:
        for bpe in calls_list:
            if curr.fields[6] == bpe.fields[6]:
                continue
            #        if bpe.c1 > curr.c1:
            #            break
            if curr.fields[0] != bpe.fields[0]:
                continue
            if (abs(int(curr.fields[1]) - int(bpe.fields[2])) < slop or abs(int(curr.fields[2]) - int(bpe.fields[1])) < slop):
                if (curr.fields[3] == bpe.fields[3]) and (abs(int(curr.fields[4]) - int(bpe.fields[5])) < slop or abs(int(curr.fields[5]) - int(bpe.fields[4])) < slop):
#                    print "removing " + bpe.shortStr() + " as a dupe of " + curr.shortStr()
                    calls_list.remove(bpe)
    
    tmp = open('calls_tmp', 'w')
    for call in calls_list:
        write_bed(call, tmp)
    
    return pybedtools.BedTool('calls_tmp')


def convert_bedpe_to_bed12(bedpe_file, track_name):
    bed12_file = open(bedpe_file + ".bed", 'w')
    subprocess.call("bedpeToBed12.py -i {0} -d 1000000000 -n \"{1}\"".format(bedpe_file, track_name), shell=True, stdout=bed12_file)
    
def uniqify(bedtool, output_dir, file_name):
    ufile = open(output_dir + "/" + file_name, 'w') 
    subprocess.call("sort -u " + output_dir + "/tmp." + file_name, shell=True, stdout=ufile) 
    os.remove(output_dir + "/tmp." + file_name)
    return pybedtools.BedTool(ufile.name)

def save_output(master_out_bed, calls, output_dir, file_name, sample_name, sv_type, seg_dups, cent_tel):
    track_name = sample_name + "_" + sv_type
    calls.saveas(output_dir + "/" + file_name + '.bedpe')
    convert_bedpe_to_bed12(output_dir + "/" + file_name + '.bedpe', track_name)    
    calls = merge_duplicate_breaks(calls, 5000)
    print sv_type + "\tNON_DUPLICATE\t" + str(len(calls))
    log.write(sv_type + "\tNON_DUPLICATE\t" + str(len(calls)) + "\n")
#    print "non-duplicate: " + str(len(calls))
    calls.saveas(output_dir + "/" + file_name + '_dedup.bedpe')
    convert_bedpe_to_bed12(output_dir + "/" + file_name + '_dedup.bedpe', track_name)            

    seg_dup_overlap = calls.pair_to_bed(seg_dups, f=1, type="either").cut(xrange(0,22)).saveas(output_dir + "/tmp." + file_name + "_dedup_segdup.bedpe")
    seg_dup_overlap = uniqify(seg_dup_overlap, output_dir, file_name + "_dedup_segdup.bedpe")
    convert_bedpe_to_bed12(output_dir + "/" + file_name + '_dedup_segdup.bedpe', track_name + "_IN_SEG_DUPS")    
    subprocess.call("cat " + output_dir + "/" + file_name + '_dedup_segdup.bedpe.bed', shell=True, stdout=master_out_bed)
    cent_tel_overlap = calls.pair_to_bed(cent_tel, f=1, type="either").cut(xrange(0,22)).saveas(output_dir + "/tmp." + file_name + "_dedup_cent_tel.bedpe")
    cent_tel_overlap = uniqify(cent_tel_overlap, output_dir, file_name + "_dedup_cent_tel.bedpe")
    convert_bedpe_to_bed12(output_dir + "/" + file_name + '_dedup_cent_tel.bedpe', track_name + "_IN_CENT_TEL")    
    subprocess.call("cat " + output_dir + "/" + file_name + '_dedup_cent_tel.bedpe.bed', shell=True, stdout=master_out_bed)

    if len(seg_dup_overlap) > 0:
        stringent_minus_sd = calls.pair_to_pair(seg_dup_overlap, type="notboth").saveas()
    else:
        stringent_minus_sd = calls
    if len(cent_tel_overlap) > 0:
        stringent_minus_ct = stringent_minus_sd.pair_to_pair(cent_tel_overlap, type="notboth").saveas()
    else:
        stringent_minus_ct = stringent_minus_sd
    
#    stringent_minus_ct.saveas(output_dir + "/stringent_minus_ct.bed")
#    stringent = calls - seg_dup_overlap - cent_tel_overlap

    print sv_type + "\tTOTAL_STRINGENT\t" + str(len(stringent_minus_ct))
    log.write(sv_type + "\tTOTAL_STRINGENT\t" + str(len(stringent_minus_ct)) + "\n")
#    print "total stringent: " + str(len(stringent_minus_ct))
    very_short_stringent = stringent_minus_ct.filter(bedpe_lt_length_filter, 1000).saveas()
    very_short_stringent.saveas(output_dir + "/" + file_name + "_dedup_stringent_very_short.bedpe")
    convert_bedpe_to_bed12(output_dir + "/" + file_name + '_dedup_stringent_very_short.bedpe', track_name + "_STRINGENT_LT_1KB")
    subprocess.call("cat " + output_dir + "/" + file_name + '_dedup_stringent_very_short.bedpe.bed', shell=True, stdout=master_out_bed)
#    stringent = stringent - very_short_stringent
    stringent_minus_vs = stringent_minus_ct.pair_to_pair(very_short_stringent, type="notboth").saveas()
    print sv_type + "\tTOTAL_STRINGENT_MINUS_VERY_SHORT\t" + str(len(stringent_minus_vs))
    log.write(sv_type + "\tTOTAL_STRINGENT_MINUS_VERY_SHORT\t" + str(len(stringent_minus_vs)) + "\n")
#    print "total stringent minus vs: " + str(len(stringent_minus_vs))

    short_stringent = stringent_minus_vs.filter(bedpe_lt_length_filter, 5000).saveas()
    short_stringent.saveas(output_dir + "/" + file_name + "_dedup_stringent_short.bedpe")
    convert_bedpe_to_bed12(output_dir + "/" + file_name + '_dedup_stringent_short.bedpe', track_name + "_STRINGENT_LT_5KB")
    subprocess.call("cat " + output_dir + "/" + file_name + '_dedup_stringent_short.bedpe.bed', shell=True, stdout=master_out_bed)
#    stringent = stringent - short_stringent
    stringent_minus_vss = stringent_minus_vs.pair_to_pair(short_stringent, type="notboth").saveas()

#    print "smv: " + str(len(stringent_minus_vss))

    stringent_high_score = stringent_minus_vss.filter(score_gte_filter, 99).saveas(output_dir + "/" + file_name + "_dedup_stringent_high_score.bedpe")
    convert_bedpe_to_bed12(output_dir + "/" + file_name + '_dedup_stringent_high_score.bedpe', track_name + "_STRINGENT_HIGH_SCORE")    
    subprocess.call("cat " + output_dir + "/" + file_name + '_dedup_stringent_high_score.bedpe.bed', shell=True, stdout=master_out_bed)
    
    stringent_low_score = stringent_minus_vss.filter(score_lt_filter, 99).saveas(output_dir + "/" + file_name + "_dedup_stringent_low_score.bedpe")
    convert_bedpe_to_bed12(output_dir + "/" + file_name + '_dedup_stringent_low_score.bedpe', track_name + "_STRINGENT_LOW_SCORE")    
    subprocess.call("cat " + output_dir + "/" + file_name + '_dedup_stringent_low_score.bedpe.bed', shell=True, stdout=master_out_bed)    

    print sv_type + "\tSEGMENTAL_DUPLICATION\t" + str(len(seg_dup_overlap))
    print sv_type + "\tIN_PERI_CENTROMERE_TELOMERE\t" + str(len(cent_tel_overlap))
    print sv_type + "\tSTRINGENT_VERY_SHORT\t" + str(len(very_short_stringent))
    print sv_type + "\tSTRINGENT_SHORT\t" + str(len(short_stringent))
    print sv_type + "\tTOTAL_STRINGENT\t" + str(len(stringent_minus_vss))
    print sv_type + "\tSTRINGENT_LOW_SCORE\t" + str(len(stringent_low_score))
    print sv_type + "\tSTRINGENT_HIGH_SCORE\t" + str(len(stringent_high_score))

    log.write( sv_type + "\tSEGMENTAL_DUPLICATION\t" + str(len(seg_dup_overlap)) + "\n")
    log.write( sv_type + "\tIN_PERI_CENTROMERE_TELOMERE\t" + str(len(cent_tel_overlap)) + "\n")
    log.write( sv_type + "\tSTRINGENT_VERY_SHORT\t" + str(len(very_short_stringent)) + "\n")
    log.write( sv_type + "\tSTRINGENT_SHORT\t" + str(len(short_stringent)) + "\n")
    log.write( sv_type + "\tTOTAL_STRINGENT\t" + str(len(stringent_minus_vss)) + "\n")
    log.write( sv_type + "\tSTRINGENT_LOW_SCORE\t" + str(len(stringent_low_score)) + "\n")
    log.write( sv_type + "\tSTRINGENT_HIGH_SCORE\t" + str(len(stringent_high_score)) + "\n")

#    print "SD: " + str(len(seg_dup_overlap)) + "; CT: " + str(len(cent_tel_overlap)) + "; Stringent (Very Short): " + str(len(very_short_stringent)) + "; Stringent (Short): " + str(len(short_stringent)) + "; Total stringent: " + str(len(stringent_minus_vss)) + "; Stringent Low Score: " + str(len(stringent_low_score)) + "; Stringent High Score: " + str(len(stringent_high_score))

hydra_calls = pybedtools.BedTool(hydra_file)
tes = pybedtools.BedTool(te_file)
common_deletions = pybedtools.BedTool(common_deletions_file)
seg_dups = pybedtools.BedTool(seg_dups_file)
cent_tel = pybedtools.BedTool(cent_tel_file)

master_out_bed = open(output_dir + "/" + sample_name + "_svs.bed", 'a')

num_calls = len(hydra_calls)
print "TOTAL\tALL\t" + str(num_calls)
log.write("TOTAL\tALL\t" + str(num_calls) + "\n")
inter_calls = hydra_calls.filter(inter_chr_filter).saveas()
print "TRANSLOCATIONS\tALL\t" + str(len(inter_calls))
log.write("TRANSLOCATIONS\tALL\t" + str(len(inter_calls)) + "\n")
save_output(master_out_bed, inter_calls, output_dir, "translocations", sample_name, "TRANSLOCATIONS", seg_dups, cent_tel)

possible_te_insertions = inter_calls.pair_to_bed(tes, f=.75).saveas()
filtered_possible_te_insertions = possible_te_insertions.filter(bedpe_reciprocal_overlap_ends_filter, 0.75).saveas()
print "TRANSLOCATIONS-POSSIBLE_TE_INSERTIONS\tALL\t" + str(len(filtered_possible_te_insertions))
log.write("TRANSLOCATIONS-POSSIBLE_TE_INSERTIONS\tALL\t" + str(len(filtered_possible_te_insertions)) + "\n")
#filtered_possible_te_insertions.saveas(output_dir + "/" + 'translocations_possible_te_insertions.bedpe')
#convert_bedpe_to_bed12(output_dir + "/" + 'translocations_possible_te_insertions.bedpe', "TRANSLOCATIONS/POSSIBLE TE INSERTIONS")
save_output(master_out_bed, filtered_possible_te_insertions, output_dir, "translocations_possible_te_insertions", sample_name, "TRANSLOCATIONS-POSSIBLE_TE_INSERTIONS", seg_dups, cent_tel)

intra_calls = hydra_calls.filter(intra_chr_filter).saveas()
#print "intra chromosomal: " + str(len(intra_calls))

expected_orientation = intra_calls.filter(expected_orientation_filter, matches=True).saveas()
#print "\texpected_orientation: " + str(len(expected_orientation))

long_indel_intra_calls = expected_orientation.filter(bedpe_gt_length_filter, insert_size).saveas()
print "DELETIONS\tALL\t" + str(len(long_indel_intra_calls))
log.write("DELETIONS\tALL\t" + str(len(long_indel_intra_calls)) + "\n")
save_output(master_out_bed, long_indel_intra_calls, output_dir, "deletions", sample_name, "DELETIONS", seg_dups, cent_tel)

possible_te_reference_insertions = long_indel_intra_calls.pair_to_bed(tes, type="ispan", f=.75).saveas()
filtered_possible_te_reference_insertions = possible_te_reference_insertions.filter(bedpe_reciprocal_overlap_ispan_filter, 0.75).saveas()
print "POSSIBLE_TE_INSERTIONS_IN_REFERENCE\tALL\t" + str(len(filtered_possible_te_reference_insertions))
log.write("POSSIBLE_TE_INSERTIONS_IN_REFERENCE\tALL\t" + str(len(filtered_possible_te_reference_insertions)) + "\n")
save_output(master_out_bed, filtered_possible_te_reference_insertions, output_dir, "possible_te_reference_insertions", sample_name, "POSSIBLE_TE_INSERTIONS_IN_REFERENCE", seg_dups, cent_tel)

common_deletions = long_indel_intra_calls.pair_to_bed(common_deletions, type="ispan", f=common_deletion_overlap_pct).saveas()
filtered_possible_common_deletions = common_deletions.filter(bedpe_reciprocal_overlap_ispan_filter, common_deletion_overlap_pct).saveas()
print "COMMON_DELETIONS\tALL\t" + str(len(filtered_possible_common_deletions))
log.write("COMMON_DELETIONS\tALL\t" + str(len(filtered_possible_common_deletions)) + "\n")
save_output(master_out_bed, filtered_possible_common_deletions, output_dir, "possible_common_deletions", sample_name, "COMMON_DELETIONS", seg_dups, cent_tel)

short_indel_intra_calls = expected_orientation.filter(bedpe_lt_length_filter, insert_size).saveas()
print "INSERTIONS\tALL\t" + str(len(short_indel_intra_calls))
log.write("INSERTIONS\tALL\t" + str(len(short_indel_intra_calls)) + "\n")
save_output(master_out_bed, short_indel_intra_calls, output_dir, "insertions", sample_name, "INSERTIONS", seg_dups, cent_tel)

unexpected_orientation = intra_calls.filter(expected_orientation_filter, matches=False).saveas()
print "INVERSION\tALL\t" + str(len(unexpected_orientation))
log.write("INVERSION\tALL\t" + str(len(unexpected_orientation)) + "\n")
save_output(master_out_bed, unexpected_orientation, output_dir, "inversions", sample_name,  "INVERSIONS", seg_dups, cent_tel)



pybedtools.cleanup()
log.close()
