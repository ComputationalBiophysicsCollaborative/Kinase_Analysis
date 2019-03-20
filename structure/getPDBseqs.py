#!/usr/bin/env python2
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
extracode = {'ABA': 'GLU', 'B3L': 'LEU', 'CAS': 'CYS', 'CIR': 'ARG', 'CME': 'CYS', 'CMT': 'CYS', 'CSD': 'CYS', 'CSO': 'CYS', 'CSS': 'CYS', 'CY0': 'CYS', 'FCL': 'PHE', 'KCX': 'LYS', 'L3O': 'LEU', 'LGY': 'LYS', 'MEA': 'PHE', 'MHO': 'MET', 'MSE': 'MET', 'NMM': 'ARG', 'OCS': 'CYS', 'OCY': 'CYS', 'PFF': 'TYR', 'PTR': 'TYR', 'SCS': 'CYS', 'SEP': 'SER', 'TPO': 'THR', 'TYI': 'TYR'}
code.update(dict((k, code[v]) for k,v in extracode.iteritems()))

# empty the output files
with open('pdbseqs.fa', 'wt') as f: 
    pass
with open('pdbseqResinds', 'wt') as f:
    pass
with open('pdbseqIndices', 'wt') as f:
    pass

for pdbfn in glob.glob(sys.argv[1] + '/*'):
    pdbname = os.path.splitext(os.path.basename(pdbfn))[0]
    print >>sys.stderr, pdbname
    
    #read pdb file raw text
    with open(pdbfn) as f:
       lines = [l.split() for l in f.readlines()] 
    seqres = [l for l in lines if l[0] == 'SEQRES'] 

    #update residue code usign MODRES
    modres = [l for l in lines if l[0] == 'MODRES'] 
    pdbcode = code.copy()
    pdbcode.update(dict((d[2],code[d[5]]) for d in modres))
    
    #also load pdb file as pdb structure
    pdb = pdbParse.loadPDB(pdbfn)

    #split up by chain
    chains = set(l[2] for l in seqres)
    for chain in chains:
        #parse SEQRES record
        cseq = [l for l in seqres if l[2] == chain]
        cseq.sort(key=lambda x: int(x[1])) #sort by SEQRES index
        rseq = [res for line in cseq for res in line[4:]]
        if rseq[0] in "AGCT": #skip DNA sequences
            print "Skip: PDB {}_{} is DNA".format(pdbname, chain)
            continue

        #trim trailing or leading weird stuff
        weird = ['ACE', 'UNK', 'NH2', 'TPO', 'AME', '69P', 'OCE', 'DAR']
        while rseq and rseq[0] in weird:
            rseq.pop(0)
        while rseq and rseq[-1] in weird:
            rseq.pop()
        #skip if there are unknown residues in the middle of the protein
        if [r for r in rseq if r == 'UNK'] != []:
            print "Skip: PDB {}_{} has unknown residues.".format(pdbname, chain)
            continue
        
        #get our final seequence
        seq = "".join(pdbcode.get(r, 'X') for r in rseq)

        #skip if sequence is short
        if len(seq) < minL: 
            print "Skip: {}_{} L = {} is < 100 for chain {}.".format(pdbname, chain, len(seq), chain)
            continue

        #now figure out map from alignment position to pdb residue id 
        #(for use in structural calculations)

        #get structure sequence from CA atoms
        c = pdb[(pdb.chain == chain)]
        CA = c[c.atom == ' CA ']
        
        #remove unknown residues and remove duplicate conformers
        fids = CA['resID'] #work with full resID to account for inserts
        filt = (CA.resName != 'UNK') & (r_[True, fids[1:] != fids[:-1]])

        #also note that there is currently a bug in numpy which causes CA.resID
        #to return a chararray, which removes whitespace. So we must access
        #resid using CA['resID'] in order to get an ndarray, which preserves
        #space. Annoying.
        
        #split up into contiguous chunks
        splitinds = where(diff(CA.resNum[filt]) > 1)[0]+1
        chunkrinds = split(CA['resID'][filt], splitinds)

        #get residue letters and split into chunks too
        CAres = array([pdbcode[r] for r in CA.resName[filt]])
        chunkres = ["".join(s) for s in split(CAres, splitinds)]

        #also keep track of 'seen' residue index (counting from 0)
        chunkinds = split(arange(len(CAres)), splitinds)
        
        #align each chunk in order
        pos = 0
        resinds, alninds = {}, {}
        for rinds, inds, res in zip(chunkrinds, chunkinds, chunkres):
            m = re.search(res, seq[pos:])
            if not m:
                raise Exception("Could not align {} inside {}, chain {}".format(res, seq, chain))
            startpos = pos + m.start()
            resinds.update(dict((startpos + n, i) for n,i in enumerate(rinds)))
            alninds.update(dict((startpos + n, i) for n,i in enumerate(inds)))
            pos += m.end()
            
        #print out full sequence
        with open('pdbseqs.fa', 'at') as f:
            print >>f, ">", pdbname+'_'+chain
            print >>f, seq
        
        #print out residue id map
        with open('pdbseqResinds', 'at') as f:
            print >>f, pdbname+'_'+chain + "|" + "|".join(str(resinds.get(i, '  -  ')) for i in range(len(seq)))
        #print out residue index map
        with open('pdbseqIndices', 'at') as f:
            print >>f, pdbname+'_'+chain,  " ".join(str(alninds.get(i, -1)) for i in range(len(seq)))
