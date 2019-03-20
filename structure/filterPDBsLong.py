#!/usr/bin/env python2
from scipy import *
import pylab
import numpy as np
import sys, os, re, errno
import pdbParse
from scipy.spatial.distance import cdist

with open(sys.argv[1]) as f:
    inlist = [l.strip() for l in f.readlines()]
with open(sys.argv[2]) as f:
    outlist = [l.strip() for l in f.readlines()]

#load map from pdb to aln index
with open(sys.argv[3]) as f:  #PDBseqResinds
    dat = [l.split('|') for l in f.readlines()]
    ridmaps = dict([(d[0],d[1:]) for d in dat])

with open(sys.argv[4]) as f: #pdbseqs_full_ID
    fullalnseqs = dict([l.split() for l in f.readlines()])

L = 175 #XXX generalize
dist = lambda x, y: sqrt(sum((x-y)**2))

def getDistances(struct):
    pdbname = struct[:4]
    chain = struct[-1]

    pdbseq2resID = ridmaps[struct] #map from pdb seq index to residueID (ie, from seqres index, to resolved residue)
    fullalnseq = fullalnseqs[struct] #seq after alignment
    
    aln2fullalnMap = where([c.isupper() or c == '-' for c in fullalnseq])[0]
    pdbseq2fullalnMap = where([c not in ".-" for c in fullalnseq])[0]
    fullaln2pdbseqMap = -ones(len(fullalnseq), dtype=int)
    fullaln2pdbseqMap[pdbseq2fullalnMap] = arange(len(pdbseq2fullalnMap))
    aln2structMap = [pdbseq2resID[fullaln2pdbseqMap[i]] if fullaln2pdbseqMap[i] != -1 else -1 for i in aln2fullalnMap]
    
    p = pdbParse.loadPDB(os.path.join('PDB', pdbname+'.pdb'))
    c = p[p.chain == chain]
    
    #this may raise an exception if salt bridge atoms missing. 
    B3Katoms = c[c.resID == aln2structMap[23]]
    aCEatoms = c[c.resID == aln2structMap[37]]
    B3K_NZ_coord = B3Katoms[B3Katoms.atom == ' NZ '][0].coords
    aCE_CD_coord = aCEatoms[aCEatoms.atom == ' CD '][0].coords
    saltDist = dist(B3K_NZ_coord, aCE_CD_coord)
    
    actDists = []
    for a,b in [(128,148), (128, 149), (128, 150)]:
        try:
            actStartAtoms = c[c.resID == aln2structMap[a]]
            actEndAtoms = c[c.resID == aln2structMap[b]]
            actDist = min(cdist(actStartAtoms.coords, actEndAtoms.coords).flatten())
        except:
            actDist = 0
        actDists.append(actDist)
    actDists = array(actDists)
    if all(actDists == 0):
        print  "              Warning, could not determine length of activation loop in {}".format(struct)

    return saltDist, actDists

def testStruct(structlist, testsalt, name):
    sdists = []
    hlens = []

    keep = []
    chk2 = []
    salt = []
    sh2 = []
    missing = []
    for struct in structlist:
        if struct not in fullalnseqs:
            print name, struct, "----   Sequence had too many gaps"
            continue

        #fullseq = fullalnseqs[struct]
        #saltDist, actDists = getDistances(struct)
        try:
            fullseq = fullalnseqs[struct]
            saltDist, actDists = getDistances(struct)
        except Exception as e:
            print name, struct, "----   Processing error (missing atom?): {}".format(repr(e))
            missing.append(struct)
            continue

        #find number of leading misaligned characters - presence of many indicates sh2 domain
        seq = re.sub('\.', '', fullseq)
        hlen = len(re.match('^([^A-Z]*)', seq).groups()[0])
        hlens.append(hlen)
        if hlen > 100:
            print name, struct, "----   Seems to have SH2 domain."
            sh2.append(struct)
            continue
        sdists.append(saltDist)

        if any(actDists > 20): #activation loop extension
            print name, struct, "----   CHK2-like crystallization: ", actDists
            chk2.append(struct)
            continue
        if testsalt and saltDist > 5.0: #B3-K to aC-glu connection
            print name, struct, "----   Matches SRC-like inactive state. Salt-bridge: {:.2f}A".format(saltDist)
            salt.append(struct)
            continue

        print name, struct
        keep.append(struct)
    return keep, chk2, salt, sh2, missing, sdists, hlens

in_keep, in_chk2, in_salt, in_sh2, in_miss, sdists, in_hlens = testStruct(inlist, True, 'In')
out_keep, out_chk2, dummy, out_sh2, out_miss, dummy, out_hlens = testStruct(outlist, False, 'Out')


with open('inList_filt', 'wt') as f:
    f.write("\n".join(sorted(in_keep)))
with open('outList_filt', 'wt') as f:
    f.write("\n".join(sorted(out_keep)))

print "Rejected PDBs:"
print ""
print "PDBS with CHK2-like conformation:"
print "In", " ".join(in_chk2)
print "Out", " ".join(out_chk2)
print ""
print "PDBS with probable SH2 domain:"
print "In", " ".join(in_sh2)
print "Out", " ".join(out_sh2)
print ""
print "In PDBS with probable salt bridge"
print " ".join(in_salt)
print ""
print "PDBS with missing atoms (probably unresolved):"
print "In", " ".join(in_miss)
print "Out", " ".join(out_miss)

sdists = array(sdists)
sdists = sdists[sdists < 99]
pylab.hist(sdists, bins=80)
pylab.xlabel('Salt Bridge Distance')

pylab.figure()
pylab.hist(in_hlens, bins=arange(0,500,2), alpha=0.5, ec='none', normed=True)
pylab.hist(out_hlens, bins=arange(0,500,2), alpha=0.5, ec='none', normed=True)
pylab.xlabel('Size of initial (SH1/2) domains')
pylab.legend(['In', 'Out'])
pylab.show()



