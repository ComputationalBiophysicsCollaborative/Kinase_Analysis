#!/usr/bin/env python
from scipy import *
import pylab
import numpy as np
import sys, os, re, errno
import pdbParse
from scipy.spatial.distance import cdist

with open(sys.argv[1]) as f:
    pdblist = [l.strip() for l in f.readlines()]

#load map from pdb to aln index
with open(sys.argv[2], 'rb') as f:  # pdbseqResinds
    dat = [l.split(b'|') for l in f.readlines()]
    ridmaps = dict([(d[0],d[1:]) for d in dat])

with open(sys.argv[3]) as f: # pdbseqs_full_ID
    fullalnseqs = dict([l.split() for l in f.readlines()])

L = 175 #XXX generalize
dist = lambda x, y: sqrt(sum((x-y)**2))

def getDistances(struct):
    pdbname = struct[:4]
    chain = struct[-1].encode('ascii')
    
    #map from pdb seq index to residueID (ie, seqres index, to resolved residue)
    pdbseq2resID = ridmaps[struct.encode('ascii')] 
    fullalnseq = fullalnseqs[struct] #seq after alignment
    
    aln2fullalnMap = where([c.isupper() or c == '-' for c in fullalnseq])[0]
    pdbseq2fullalnMap = where([c not in ".-" for c in fullalnseq])[0]
    fullaln2pdbseqMap = -ones(len(fullalnseq), dtype=int)
    fullaln2pdbseqMap[pdbseq2fullalnMap] = arange(len(pdbseq2fullalnMap))
    aln2structMap = [pdbseq2resID[fullaln2pdbseqMap[i]] 
                  if fullaln2pdbseqMap[i] != -1 else -1 for i in aln2fullalnMap]
    
    p = pdbParse.loadPDB(os.path.join('PDB', pdbname+'.pdb'))
    c = p[p.chain == chain]
    
    actDists = []
    for a,b in [(128,148), (128, 149), (128, 150)]:
        try:
            actStartAtoms = c[c.resID == aln2structMap[a]]
            actEndAtoms = c[c.resID == aln2structMap[b]]
            actDist = min(cdist(actStartAtoms.coords, actEndAtoms.coords
                         ).flatten())
        except:
            actDist = 0
        actDists.append(actDist)
    actDists = array(actDists)
    if all(actDists == 0):
        print("              Warning, could not determine length of "
              "activation loop in {}".format(struct))

    return actDists

def testStruct(structlist, name):
    sdists = []
    hlens = []

    keep = []
    chk2 = []
    sh2 = []
    missing = []
    for struct in structlist:
        if struct not in fullalnseqs:
            print(name, struct, "----   Sequence had too many gaps")
            continue

        fullseq = fullalnseqs[struct]
        actDists = getDistances(struct)

        #find number of leading misaligned characters - presence of many
        #indicates sh2 domain
        seq = re.sub('\.', '', fullseq)
        hlen = len(re.match('^([^A-Z]*)', seq).groups()[0])
        hlens.append(hlen)
        if hlen > 100:
            print(name, struct, "----   Seems to have SH2 domain.")
            sh2.append(struct)
            continue

        if any(actDists > 20): #activation loop extension
            print(name, struct, "----   CHK2-like crystallization: ",
                 "({})".format(", ".join("{:.2f}".format(x) for x in actDists)))
            chk2.append(struct)
            continue

        print(name, struct)
        keep.append(struct)
    return keep, chk2, sh2, missing, sdists, hlens

keep, chk2, sh2, miss, dummy, hlens = testStruct(pdblist, 'PDB')

with open('pdblist_filt', 'wt') as f:
    f.write("\n".join(keep))

print("Rejected PDBs:")
print("")
print("PDBS with CHK2-like conformation:")
print(" ".join(chk2))
print("")
print("PDBS with probable SH2 domain:")
print(" ".join(sh2))
print("")
print("PDBS with missing atoms (probably unresolved):")
print(" ".join(miss))

pylab.figure()
pylab.hist(hlens, bins=arange(0,500,2), alpha=0.5, ec='none', normed=True)
pylab.xlabel('Size of initial (SH1/2) domains')
pylab.show()



