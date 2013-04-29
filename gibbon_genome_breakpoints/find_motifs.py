from Bio import SeqIO
import motility
import sys
import re

ref = sys.argv[1]
motif_list = sys.argv[2]

sys.stderr.write("loading sequences\n")
record_dict = SeqIO.index(ref,"fasta")

sys.stderr.write("loaded " + str(len(record_dict.keys())) + " sequences\n") 

def re_from_iupac(motif):
    iupacdict = {'A':'A',
                 'C':'C',
                 'G':'G',
                 'T':'T',
                 'M':'[AC]',
                 'R':'[AG]',
                 'W':'[AT]',
                 'S':'[CG]',
                 'Y':'[CT]',
                 'K':'[GT]',
                 'V':'[ACG]',
                 'H':'[ACT]',
                 'D':'[AGT]',
                 'B':'[CGT]',
                 'X':'[ACGT]',
                 'N':'[ACGT]'}
    regex = ''
    for c in motif:
        regex += iupacdict[c]
    return regex

motifs = {}
for line in open(motif_list).readlines():
    name, motif = line.rstrip().split("\t")
    motifs[name] = motif

for motif_name in motifs.keys():
    out_file = open(motif_name + ".bed", "w")
    motif = motifs.get(motif_name)
    
    for chrom_name in record_dict.keys():
        sys.stderr.write("searching for motif " + motif + " in chrom " + chrom_name +"\n")
        myseq = record_dict[chrom_name].seq
        seqlen = len(str(myseq))
        sys.stderr.write("seqlen is " + str(seqlen) + "\n")
        motif_re = re.compile(re_from_iupac(motif))
        sys.stderr.write("searching forward\n")
        for i in motif_re.finditer(str(myseq)):
            out_file.write(chrom_name + "\t" + str(i.start()) + "\t" + str(i.end()) + "\t" + str(myseq)[i.start():i.end()] + "\t1\t+\n")
        sys.stderr.write("searching rev comp\n")
        for i in motif_re.finditer(str(myseq.reverse_complement())):
            out_file.write(chrom_name + "\t" + str(seqlen - i.end()) + "\t" + str(seqlen - i.start()) + "\t" + str(myseq)[(seqlen - i.end()):(seqlen - i.start())] + "\t1\t-\n")
    out_file.close()
