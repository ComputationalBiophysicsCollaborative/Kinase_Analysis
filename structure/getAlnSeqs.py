#!/usr/bin/env python
from scipy import *
import sys, re, argparse

parser = argparse.ArgumentParser(description='Extract Seqs')
parser.add_argument('stkfile')
parser.add_argument('mode', choices=['orig', 'full', 'aln', 'all', 'pdb'])
parser.add_argument('--showid', action='store_true')
parser.add_argument('--Xgap', action='store_true')
parser.add_argument('--unalign') #make these positions unaligned even if they align by default
args = parser.parse_args(sys.argv[1:])

pdbseqs = {}
pdbfns = []
with open(args.stkfile) as f:
    dat = [l.split() for l in f.readlines()[:-1]
           if (not l.isspace()) and (not l.startswith('#'))]
    for x in dat:
        if x[0] in pdbfns:
            pdbseqs[x[0]] += x[1]
        else:
            pdbfns.append(x[0])
            pdbseqs[x[0]] = x[1]

mode = args.mode
if args.unalign:
    parts = [l.split('-') for l in args.unalign.split(',')]
    unalign = concatenate([range(int(p[0]),int(p[1])+1)
                           if len(p) > 1 else [int(p)] for p in parts])

for pdbfn in pdbfns:
    seq = pdbseqs[pdbfn]
    head = pdbfn + ' ' if args.showid else ''

    if args.Xgap:
        seq = seq.replace('X', '-')

    if args.unalign:
        matches = [m.start() for m in re.finditer('[A-Z-]', seq)]
        ls = list(seq)
        for i in unalign:
            c = ls[matches[i]]
            ls[matches[i]] = c.lower() if c != '-' else '.'
        seq = "".join(ls)

    if mode == 'orig':
        print(head + seq)
    elif mode == 'full':
        print(head + re.sub("[xX.-]", "", seq))
    elif mode == 'aln':
        print(head + re.sub("[a-z.]", "", seq))
    elif mode == 'pdb':
        print(head + re.sub("[.-]", "", seq))
    elif mode == 'all':
        print(head + seq)

