#!/usr/bin/env python
from math import sqrt as msqrt
from scipy import *
import pdbParse
import sys, os, glob, errno
from scipy.spatial.distance import cdist
import time

code = pdbParse.residueCode
extracode = {'ABA': 'GLU', 'B3L': 'LEU', 'CAS': 'CYS', 'CIR': 'ARG', 'CME':
             'CYS', 'CMT': 'CYS', 'CSD': 'CYS', 'CSO': 'CYS', 'CSS': 'CYS',
             'CY0': 'CYS', 'FCL': 'PHE', 'KCX': 'LYS', 'L3O': 'LEU', 'LGY':
             'LYS', 'MEA': 'PHE', 'MHO': 'MET', 'MSE': 'MET', 'NMM': 'ARG',
             'OCS': 'CYS', 'OCY': 'CYS', 'PFF': 'TYR', 'PTR': 'TYR', 'SCS':
             'CYS', 'SEP': 'SER', 'TPO': 'THR', 'TYI': 'TYR'}
code.update(dict((k, code[v]) for k,v in extracode.items()))
code = {k.encode('ascii'): v for k,v in code.items()}
#note that these extra codes are unreliable, since some PDBS assign the same
#name to different residues. So correct for that using MODRES below

PDBdir = sys.argv[1]

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

mkdir_p('distancesNHA')
mkdir_p('distancesCA')
mkdir_p('distancesCB')

def vdist(a,b):
    return sqrt(sum((a-b)**2))


for pdbfn in glob.glob(PDBdir + '/*'):
    pdb = os.path.splitext(os.path.basename(pdbfn))[0]
    print(pdb)
    p = pdbParse.loadPDB(pdbfn)
    if len(p) == 0:
        print("Empty file? ", pdbfn)
    chains = [c for c in set(p.chain)]

    for chain in chains:
        name = pdb+'_'+chain.decode('ascii')

        if os.path.exists(os.path.join('distances', name+'.npy')):
            print("Already done {}".format(name))
            continue

        #get atoms, but remove unknown residues and remove duplicate conformers
        #Important that this is done in exactly the same way as in
        #getPDBseq.py so that "seen" residue indexes match
        c = p[(p.chain == chain)]
        c = c[((c.altLoc == b' ') | (c.altLoc == b'A')) & (c.resName != b'UNK')]

        CA = c[c.atom == b' CA ']

        resids = CA.resID
        L = len(resids)
        pairs = [(i,j) for i in range(L-1) for j in range(i+1,L)]


        # compute Nearest-heavy-atom (NHA) distances
        atoms = c[c.element != b' H']
        atomids = atoms.resID
        cc = [atoms.coords[argwhere(atomids == id).ravel()] for id in resids]
        NHAdist = array([min(cdist(cc[i], cc[j]).ravel()) for i,j in pairs])

        # compute CA distance
        CAdist = cdist(CA.coords, CA.coords)[triu_indices(L,k=1)]

        # compute CB distance.  Count CA as CB for gly
        CB = c[(c.atom == b' CB ') |
               ((c.atom == b' CA ') & (c.resName == b'GLY'))]
        inds = [argwhere(CB.resID == id).ravel() for id in resids]
        cc = array([CB.coords[ind[0]] if len(ind) > 0 else [nan,nan,nan]
                    for ind in inds])
        CBdist = cdist(cc, cc)[triu_indices(L,k=1)]

        # sanity checks
        assert(len(CAdist) == len(CBdist))
        assert(len(CAdist) == len(NHAdist))
        # sanity check for clashes
        clashes = [(resids[i], resids[j]) for n,(i,j) in enumerate(pairs)
                   if abs(i-j) > 4 and NHAdist[n] < 2.0]
        # 2.0 isdistance at which pymol draws bonds
        if len(clashes) > 0:
            clashes = [(x.decode('ascii'),y.decode('ascii')) for x,y in clashes]
            print('Clashes for pdb {} in residues {}'.format(name, clashes))

        save(os.path.join('distancesNHA', name), NHAdist)
        save(os.path.join('distancesCA', name), CAdist)
        save(os.path.join('distancesCB', name), CBdist)
