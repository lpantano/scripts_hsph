#script to detect modificatino in read2 of tail-seq
import gzip
import re
import numpy
import os
from collections import Counter
from utils import open_fastq, file_transaction

w = {}
w['A'] = -8
w['C'] = -8
w['G'] = -8
w['T'] = 2
w['N'] = -5

def getsc(nt):
        if w.has_key(nt):
                return(w[nt])
        else:
                return -10


def poly_A_percentage(seq):
    predicted = {}
    if re.search("TTT", seq):
        for s in range(30):
            #print s
            if seq[s:].startswith("TT"):
            #string[s]=="T":
                for e in range(s+5, s+60, 2):
                    #print e
                    sc = sum([1 for nt in seq[s:e] if nt == "T"])
                    scoreT = 1.0 * sc / (e-s)
                    #print "%s %s %s" % (s, e, scoreT)
                    if scoreT > 0.70:
                        predicted[(e-s+1, scoreT)] = (s, e)
    if len(predicted.keys()) > 0:
        sorted_scores = sorted(predicted.keys(), reverse=True)
        #print "max score %s pos %s " % (sorted_scores[0], predicted[sorted_scores[0]])
        return predicted[sorted_scores[0]]
    else:
        False


def polyA(string):
    smax = -1
    emax = -1
    scmax = -100
    if re.search("TTT", string):
        for s in range(30):
            sc = 0
            #print s
            if re.match("TT", string[s:]):
            #string[s]=="T":
                for e in range(s+5, 50-s, 2):
                    #print e
                    sc = numpy.sum(map(getsc, string[s:(e)]))
                    nT = string[s:(e)].count("T")
                    nO = len(string[s:(e)])-nT
                    if nO > 0:
                #       print "other N %s" % nO
                        sc = nT - nO + sc
                    #print string[s:(s+e)] 
                    #print sc
                    #print "%s %s %s %s %s" % (sc,scmax,s,e,string[s:(e)])
                    if string[e] == "T" and string[e-1] != "T":
                        sc += 2
                        if sc >= scmax+6:
                            scmax = sc
                            smax = s
                            emax = e+1
                            if s > 0:
                                if string[s-1] == "T":
                                    smax = s-1
                                    scmax += 2
                    if string[e-2] == "T" and string[e-1] != "T":
                        sc += 8
                        if sc >= scmax+6:
                            scmax = sc
                            smax = s
                            emax = e-1
                            if s > 0:
                                if string[s-1] == "T":
                                    smax = s-1
                                    scmax += 2
                    if sc >= scmax + 6 and string[e-1] == "T":
                        scmax = sc
                        smax = s
                        emax = e
                        if s > 0:
                            if string[s-1] == "T":
                                smax = s-1
                                scmax += 2
                    #print "# %s %s %s" % (scmax,smax,emax)
        if scmax > 0:
            return([smax, emax])
    else:
            return(False)


def _adapter(seq, qual):
    """detect 5nt TAG (GTCAG) and remove from sequence.
    It should be aroung 15-22 position in read
    """
    TAG = "GTCAG"
    tag_pos = seq.find(TAG, 15, 25)
    #print tag_pos
    if tag_pos:
        tag_pos += 5
        return seq[tag_pos:], qual[tag_pos:]
    else:
        return False


def _test():
        s1 = "CCCCGCATTAAACTTGTCAGAACCAGAGTNATCTTTTTTTTATTTTTTATCTTTTTGATTTATTTTCAGCTCTTCTTTTTCAGTCAAGAATTCTTGCTATTAGGAAAATAATTCCAGATACCATTATAGTAAATATTGCTAAAATGCAAAATACTAATAAAACCTTAGTAAAGTATGAAACTAAAACTAATAGGAAAATTAGAATTGGTGATAATGCTGATAATGAACAATATGAAGTAA"
        seq, qual = _adapter(s1, s1)
        print seq
        ns = poly_A_percentage(seq)
        print ns
        print seq[:ns[0]]
        print seq[ns[0]:ns[1]]


def detect(in_file, out_prefix):
    out_name = out_prefix + "-polyA.dat.gz"
    out_name_false = out_prefix + "-none.dat.gz"
    counts = Counter()
    print "reading file %s" % in_file
    print "creating files %s %s" % (out_name, out_name)
    if os.path.exists(out_name):
        return out_name
    with file_transaction(out_name) as tx_out_file:
        with open_fastq(in_file) as handle, gzip.open(tx_out_file, 'w') as out, gzip.open(out_name_false, 'w') as out_false:
            for line in handle:
                #print line
                if line.startswith("@HISEQ"):
                    #print line
                    name = line.strip()
                    seq = handle.next().strip()
                    handle.next().strip()
                    qual = handle.next().strip()
                    find = _adapter(seq, qual)
                    #print "%s %s" % (seq, find)
                    if find:
                        seq, qual = find
                        ns = poly_A_percentage(seq)
                        if ns:
                            if ns[1]-ns[0] >= 6:
                                #print "positions are" + str(ns[0]) + ".." + str(ns[1])
                                mod = seq[:ns[0]]
                                seq_polyA = seq[ns[0]:ns[1]]
                                seq_gene = seq[ns[1]:]
                                qual_polyA = qual[ns[0]:ns[1]]
                                qual_gene = qual[ns[1]:]
                                #print "%s\t%s\t%s\t%s\t%s\t%s\n" % (name,mod,sf,qf)
                                out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (name, ns[0], ns[1], mod, seq_polyA, qual_polyA, seq_gene, qual_gene))
                                counts['polyA'] += 1
                                if len(mod) > 0:
                                    counts['mod'] += 1
                            else:
                                counts['shortA'] += 1
                                out_false.write("%s\t%s\t%s\t%s\n" % ("shortA", name, seq, qual))
                        else:
                            counts['noA'] += 1
                            out_false.write("%s\t%s\t%s\t%s\n" % ("None", name, seq, qual))
                    else:
                        out_false.write("%s\t%s\t%s\t%s\n" % ("No_tag", name, seq, qual))
                        counts['notag'] += 1
        with open(out_name + ".stat", 'w') as handle:
            handle.write("%s" % counts)
    return out_name
