#!/usr/bin/env python2
from scipy import *
import sys

with open(sys.argv[1]) as f: #pdbseqResinds
    dat = [l.split('|') for l in f.readlines()]
    inds = dict([(d[0], d[1:]) for d in dat])

with open(sys.argv[2]) as f: #pdbseqs_full_ID
    alnseqs = dict([l.split() for l in f.readlines()])

for pdb in alnseqs:
    pdbpos = where([c not in '.-' for c in alnseqs[pdb]])[0]
    indmap = dict(zip(pdbpos, inds[pdb]))
    alnpos = where([c.isupper() or c == '-' for c in alnseqs[pdb]])[0]
    print pdb + "|" + '|'.join(indmap.get(i, '  -  ') for i in alnpos)
