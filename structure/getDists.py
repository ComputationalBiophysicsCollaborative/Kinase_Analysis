#!/usr/bin/env python2
from math import sqrt as msqrt
from scipy import *
import pdbParse
import sys, os, glob, errno
from scipy.spatial.distance import cdist

code = pdbParse.residueCode
extracode = {'ABA': 'GLU', 'B3L': 'LEU', 'CAS': 'CYS', 'CIR': 'ARG', 'CME': 'CYS', 'CMT': 'CYS', 'CSD': 'CYS', 'CSO': 'CYS', 'CSS': 'CYS', 'CY0': 'CYS', 'FCL': 'PHE', 'KCX': 'LYS', 'L3O': 'LEU', 'LGY': 'LYS', 'MEA': 'PHE', 'MHO': 'MET', 'MSE': 'MET', 'NMM': 'ARG', 'OCS': 'CYS', 'OCY': 'CYS', 'PFF': 'TYR', 'PTR': 'TYR', 'SCS': 'CYS', 'SEP': 'SER', 'TPO': 'THR', 'TYI': 'TYR'}
code.update(dict((k, code[v]) for k,v in extracode.iteritems()))
#note that these extra codes are unreliable, since some PDBS assign the same 
#name to different residues. So correct for that using MODRES below

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

mkdir_p('distances')
mkdir_p('distancesCA')
mkdir_p('distancesCB')

def vdist(a,b):
    return sqrt(sum((a-b)**2))

for pdbfn in glob.glob(sys.argv[1] + '/*'):
#for pdbfn in [sys.argv[1] + '/' + x + '.pdb' for x in ['1IEP', '2GQG']]:
    pdb = os.path.splitext(os.path.basename(pdbfn))[0]
    print >>sys.stderr, pdb
    p = pdbParse.loadPDB(pdbfn)
    chains = list(set(p.chain))
    for chain in chains:
        if os.path.exists(os.path.join('distances', pdb+'_'+chain+'.npy')):
            print "Already done {}".format(pdb+'_'+chain)
            continue

        c = p[(p.chain == chain)]
        CA = c[c.atom == ' CA ']
        CB = c[(c.atom == ' CB ') | ((c.atom == ' CA ') & (c.resName == 'GLY'))] #for CB calculation, count CA as CB for gly

        #remove unknown residues and remove duplicate conformers
        #Important that this is done in exactly the same way as in 
        #getPDBseq.py so that "seen" residue indexes match
        fidsA = CA['resID'] #work with full resID to account for inserts
        filtA = (CA.resName != 'UNK') & (r_[True, fidsA[1:] != fidsA[:-1]])

        fidsB = CB['resID'] #work with full resID to account for inserts
        filtB = (CB.resName != 'UNK') & (r_[True, fidsB[1:] != fidsB[:-1]])

        resids = CA['resID'][filtA]
        nres = len(resids)
        if nres == 0:
            continue

        def findfirst(cond):
            res = argwhere(cond).flatten()
            return res[0] if len(res) > 0 else None

        def findall(cond):
            return argwhere(cond).flatten()

        if 1:#not os.path.exists(os.path.join('distances', pdb+'_'+chain+'.npy')):
            rinds = [findall((c['resID'] == resids[i]) & (c.atom.element != ' H')) for i in range(nres)]
            residues = [c[i].coords for i in rinds]
            distances = array([min(cdist(residues[i], residues[j]).flatten()) for i in range(nres-1) for j in range(i+1, nres)])

            rinds = [findfirst((CA['resID'] == resids[i]) & filtA) for i in range(nres)]
            residues = [CA[i].coords for i in rinds]
            CAdist = array([vdist(residues[i], residues[j]) for i in range(nres-1) for j in range(i+1, nres)])
            
            rinds = [findfirst((CB['resID'] == resids[i]) & filtB) for i in range(nres)]
            residues = [CB[i].coords if i != None else nan for i in rinds]
            CBdist = array([vdist(residues[i], residues[j]) for i in range(nres-1) for j in range(i+1, nres)])

            assert(len(CAdist) == len(CBdist))
            assert(len(CAdist) == len(distances))

            save(os.path.join('distances', pdb+'_'+chain), distances)
            save(os.path.join('distancesCA', pdb+'_'+chain), CAdist)
            save(os.path.join('distancesCB', pdb+'_'+chain), CBdist)
