#!/usr/bin/env python

import sys
from Bio import SeqIO

def countBadBases(seq1, seq2, pos, length):
    n1 = 0
    d1 = 0
    n2 = 0
    d2 = 0
    for i in range(length):
        b = seq1[pos]
        if b in "Nn":
            n1 += 1
        elif b == '-':
            d1 += 1
        b = seq2[pos]
        if b in "Nn":
            n2 += 1
        elif b == '-':
            d2 += 1
        pos += 1
    return (n1, d1, n2, d2)

def getBounds(seq1, seq2, length=11, maxd=100, maxn=2):
    seqlen = min(len(seq1), len(seq2))
    k1 = 0
    k2 = seqlen - length
    good1 = False
    good2 = False

#    print seq1
#    print seq2
#    print len(seq1)
#    print len(seq2)
#    raw_input()
    for p in range(maxd):
        n1, d1, n2, d2 = countBadBases(seq1, seq2, p, length)
        if d1 == 0 and d2 == 0 and n1 <= maxn and n2 <= maxn:
            k1 = p
            good1 = True
            break

    for p2 in range(maxd):
        p = seqlen - length - p2 
        n1, d1, n2, d2 = countBadBases(seq1, seq2, p, length)
        if d1 == 0 and d2 == 0 and n1 <= maxn and n2 <= maxn:
            k2 = p
            good2 = True
            break

    if good1 and good2:
        return (k1, k2)
    else:
        return None

def writeMfile(mout, scffile, alnfile, imgfile, length=11, maxd=150):
    handle = SeqIO.parse(alnfile, "fasta")
    ref = next(handle).seq
    read = next(handle).seq
    bounds = getBounds(ref, read, length=length, maxd=maxd)
    if bounds:
        sys.stderr.write("Bounds found at positions {}, {}.\n".format(bounds[0], bounds[1]))
        mout.write("bisulfitePCRSeqAnalysis('{}', '{}', '{}', '{}', '{}', '{}', '{}');\n\n".format(
            scffile, ref,
            ref[bounds[0] : bounds[0] + length - 1],
            ref[bounds[1] : bounds[1] + length - 1],
            read[bounds[0] : bounds[0] + length - 1],
            read[bounds[1] : bounds[1] + length - 1],
            imgfile))

if __name__ == "__main__":
    args = sys.argv[1:]
    mfile = args[0]
    infiles = args[1:]
    nfiles = len(infiles)
    if nfiles % 3 == 0:
        with open(mfile, "w") as mout:
            for i in range(0, nfiles, 3):
                writeMfile(mout, infiles[i], infiles[i+1], infiles[i+2])
