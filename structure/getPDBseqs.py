#!/usr/bin/env python
from math import sqrt as msqrt
from scipy import *
import pdbParse
import sys, os, glob, re
from scipy.spatial.distance import cdist

minL = int(sys.argv[2])
code = pdbParse.residueCode

#some PDBs have nonstandard residues. PDBs are supposed to list them in MODRES
#records, but don't always. So we keep a "default" list of nonstandard residues
#which is overridden if MODRES is found (important since PDBs sometimes assign
#the same name to different residues)
extracode = {'ABA': 'GLU', 'B3L': 'LEU', 'CAS': 'CYS', 'CIR': 'ARG', 'CME':
             'CYS', 'CMT': 'CYS', 'CSD': 'CYS', 'CSO': 'CYS', 'CSS': 'CYS',
             'CY0': 'CYS', 'FCL': 'PHE', 'KCX': 'LYS', 'L3O': 'LEU', 'LGY':
             'LYS', 'MEA': 'PHE', 'MHO': 'MET', 'MSE': 'MET', 'NMM': 'ARG',
             'OCS': 'CYS', 'OCY': 'CYS', 'PFF': 'TYR', 'PTR': 'TYR', 'SCS':
             'CYS', 'SEP': 'SER', 'TPO': 'THR', 'TYI': 'TYR'}
code.update(dict((k, code[v]) for k,v in extracode.items()))
code = {k.encode('ascii'): v for k,v in code.items()}

# empty the output files
with open('pdbseqs.fa', 'wt') as f:
    pass
with open('pdbseqResinds', 'wt') as f:
    pass
with open('pdbseqIndices', 'wt') as f:
    pass

for pdbfn in glob.glob(sys.argv[1] + '/*'):
    pdbname = os.path.splitext(os.path.basename(pdbfn))[0]
    print(pdbname)

    #read pdb file raw text
    with open(pdbfn, 'rb') as f:
       lines = [l.split() for l in f.readlines()]
    seqres = [l for l in lines if l[0] == b'SEQRES']

    #update residue code usign MODRES
    modres = [l for l in lines if l[0] == b'MODRES']
    pdbcode = code.copy()
    pdbcode.update(dict((d[2],code[d[5]]) for d in modres))

    #also load pdb file as pdb structure
    pdb = pdbParse.loadPDB(pdbfn)

    #split up by chain
    chains = set(l[2] for l in seqres)
    for chain in chains:
        name = pdbname+'_'+chain.decode('ascii')

        #parse SEQRES record
        cseq = [l for l in seqres if l[2] == chain]
        cseq.sort(key=lambda x: int(x[1])) #sort by SEQRES index
        rseq = [res for line in cseq for res in line[4:]]
        if rseq[0] in b"AGCT": #skip DNA sequences
            print("Skip: PDB {}_{} is DNA".format(pdbname, chain))
            continue

        #trim trailing or leading weird stuff
        weird = [b'ACE', b'UNK', b'NH2', b'TPO', b'AME', b'69P', b'OCE', b'DAR']
        while rseq and rseq[0] in weird:
            rseq.pop(0)
        while rseq and rseq[-1] in weird:
            rseq.pop()
        #skip if there are unknown residues in the middle of the protein
        if [r for r in rseq if r == b'UNK'] != []:
            print("Skip: PDB {}_{} has unknown residues.".format(pdbname,chain))
            continue

        #get our final sequence
        seq = "".join(pdbcode.get(r, 'X') for r in rseq)
        L = len(seq)

        #skip if sequence is short
        if len(seq) < minL:
            print("Skip: {} L = {} is < 100 for chain {}.".format(
                  name, len(seq), chain.decode('ascii')))
            continue

        #now figure out map from alignment position to pdb residue id
        #(for use in structural calculations)

        #get structure sequence from CA atoms
        c = pdb[(pdb.chain == chain)]
        c = c[((c.altLoc == b' ') | (c.altLoc == b'A')) & (c.resName != b'UNK')]
        CA = c[c.atom == b' CA ']

        #also note that there is currently a bug in numpy which causes CA.resID
        #to return a chararray, which removes whitespace. So we must access
        #resid using CA['resID'] in order to get an ndarray, which preserves
        #space. Annoying.

        #split up into contiguous chunks
        splitinds = where(diff(CA.resNum) > 1)[0]+1
        chunkrinds = split(CA.resID, splitinds)

        #get residue letters and split into chunks too
        CAres = array([pdbcode[r] for r in CA.resName])
        chunkres = ["".join(s) for s in split(CAres, splitinds)]

        #also keep track of 'seen' residue index (counting from 0)
        chunkinds = split(arange(len(CAres)), splitinds)

        #align each chunk in order
        pos = 0
        resinds, alninds = {}, {}
        for rinds, inds, res in zip(chunkrinds, chunkinds, chunkres):
            m = re.search(res, seq[pos:])
            if not m:
                raise Exception("Could not align {} inside {}, chain {}".format(
                                                               res, seq, chain))
            startpos = pos + m.start()
            resinds.update(dict((startpos + n, i) for n,i in enumerate(rinds)))
            alninds.update(dict((startpos + n, i) for n,i in enumerate(inds)))
            pos += m.end()

        #print out full sequence
        with open('pdbseqs.fa', 'at') as f:
            print(">", name, file=f)
            print(seq, file=f)

        #print out residue id map
        with open('pdbseqResinds', 'at') as f:
            idmap = "|".join(resinds.get(i, b'  -  ').decode('ascii') for i in range(L))
            print(name + "|" + idmap, file=f)

        #print out residue index map
        with open('pdbseqIndices', 'at') as f:
            ridmap = " ".join(str(alninds.get(i, -1)) for i in range(L))
            print(name, ridmap, file=f )
