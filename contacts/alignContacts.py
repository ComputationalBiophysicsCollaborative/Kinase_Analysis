#!/usr/bin/env python2
from scipy import *
import scipy.io as sio
import matplotlib
#matplotlib.use('Agg')
import pylab
import matplotlib.patheffects as PathEffects
import matplotlib.transforms as transforms
import numpy as np
import sys, os, re, errno
from itertools import groupby
from annotate import annotate

scriptPath = os.path.dirname(os.path.realpath(__file__))

outpath = sys.argv[1]
distpath = sys.argv[2]

cutoff = 10

#load map from pdb to aln index
with open(sys.argv[3]) as f:  #PDBseqIndices
    dat = [l.split() for l in f.readlines()]
    indmaps = dict([(d[0], array([int(x) for x in d[1:]])) for d in dat])

with open(sys.argv[4]) as f: #PDBseqsHHfull_filtergap
    fullalnseqs = dict([l.split() for l in f.readlines()])

with open(sys.argv[5]) as f:
    dat = [l.split(',') for l in f.readlines()]
    regions = [(n,int(s),int(e)) for n,s,e in dat]
lobediv = [e for n,s,e in regions if n == 'Hinge'][0]

def fillGroups(positions, L, clr):
    for k, g in groupby(enumerate(positions), lambda (i,x):i-x):
        inds = [x[1] for x in g]
        start, end = inds[0], inds[-1]+1
        pylab.fill([start+0.5, start+0.5, end+0.5, end+0.5], [+0.5,L+0.5,L+0.5,+0.5],color=clr, ec='none', lw=0)
        pylab.fill([+0.5,L+0.5,L+0.5,+0.5],[start+0.5, start+0.5, end+0.5, end+0.5], color=clr, ec='none', lw=0)

def plotContacts(name):
    pdbseq2structMap = indmaps[name] #map from pdb seq index to structural index (ie, from seqres index, to resolved residue #)
    fullalnseq = fullalnseqs[name] #seq after alignment
    
    aln2fullalnMap = where([c.isupper() or c == '-' for c in fullalnseq])[0]
    pdbseq2fullalnMap = where([c not in ".-" for c in fullalnseq])[0]
    fullaln2pdbseqMap = -ones(len(fullalnseq), dtype=int)
    fullaln2pdbseqMap[pdbseq2fullalnMap] = arange(len(pdbseq2fullalnMap))
    aln2structMap = [pdbseq2structMap[fullaln2pdbseqMap[i]] if fullaln2pdbseqMap[i] != -1 else -1 for i in aln2fullalnMap]

    pdbseq = [x for x in fullalnseq if x not in ".-"]
    alnseq = [x for x in fullalnseq if x not in "." and not x.islower()]

    distances = load(os.path.join(distpath, name + '.npy'))
    L = int(((1+sqrt(1+8*len(distances)))/2) + 0.5) 

    #form matrix
    distancesM = zeros((L,L),dtype=float)
    distancesM[triu_indices(L,k=1)] = distances
    distancesM = distancesM + distancesM.T

    ########################################
    #first, plot the raw contact map for the entire sequence 
    #that is, the contact map size will be equal to the length of the seqres header (minus junk)
    #shows inserts as red, missing residues as blue, and deletes as dotted lines (or colored lines?)

    pylab.figure(figsize=(12,12))
    ax = pylab.axes() 
    trans = transforms.blended_transform_factory(ax.transAxes, ax.transData)

    sL = len(pdbseq2structMap)
    distmat = zeros((sL,sL))*nan
    for n,(i,j) in enumerate([(i,j) for i in range(sL-1) for j in range(i+1,sL)]):
        if pdbseq2structMap[i] == -1 or pdbseq2structMap[j] == -1:
            continue
        distmat[i,j] = distancesM[pdbseq2structMap[i], pdbseq2structMap[j]]
        distmat[j,i] = distancesM[pdbseq2structMap[j], pdbseq2structMap[i]]
    contacts = distmat < cutoff
    contacts[diag([i != -1 for i in pdbseq2structMap])] = True
    
    cim = ones(distmat.shape + (3,))
    cim[contacts] = ones(3)*0.4
    pylab.imshow(cim, origin='lower', extent=(+0.5,sL+0.5,+0.5,sL+0.5), interpolation='nearest')

    #annotate inserts relative to alignment (red)
    inserts = where([c.islower() for c in pdbseq])[0]
    fillGroups(inserts, sL, (0.7,0.4,0.4,0.3))

    #annotate missign residues (blue)
    missing = where([i == -1 for i in pdbseq2structMap])[0]
    fillGroups(missing, sL, (0.4,0.4,0.7,0.3))

    #map from alignment index to *closest* pdb seq index
    remap = lambda i,s: searchsorted(pdbseq2fullalnMap, aln2fullalnMap[i], side=s)

    #annotate regions
    xregions = [(n, remap(s, 'left'), remap(e-1, 'left')+1) for n,s,e in regions]
    annotate(xregions, remap(lobediv, 'left'), sL, pylab.gca(), zorder=1)

    #plot lines at deletion points (green)
    deletepos = searchsorted(pdbseq2fullalnMap, where([x == '-' for x in fullalnseq])[0])
    #for d in set(deletepos):
    #    pylab.axvline(d+0.5, color='g', zorder=9)
    #    pylab.axhline(d+0.5, color='g', zorder=9)
    transY = transforms.blended_transform_factory(ax.transAxes, ax.transData)
    altern = False
    for d, g in groupby(deletepos):
        pylab.axvline(d+0.5, color='g', zorder=9)
        pylab.axhline(d+0.5, color='g', zorder=9)
        pylab.text(0.005 + altern*0.01, d+0.5, str(len(list(g))), verticalalignment='center', horizontalalignment='right', transform=transY, path_effects=[PathEffects.withStroke(linewidth=3, foreground="w")], color='g', zorder=100)
        altern = not altern

    #plot seqeunce
    altern = 0
    for n,c in enumerate(pdbseq):
        pylab.text(1.02 + altern*0.014, n+1, c, verticalalignment='center', horizontalalignment='center', transform=transY, zorder=100)
        altern = (altern+1)%6
    
    #plot dotted line around alignment region
    start = remap(0, 'left')-0.5
    end = remap(-1, 'right')+1.5
    pylab.fill([start, start, end, end],[start, end, end, start], fill=False, ls='dashed', zorder=10)

    pylab.title(name + ", full PDB sequence ({}A)".format(cutoff))
    pylab.subplots_adjust(bottom=0.05, right=0.90, top=0.95, left=0.05)

    ########################################
    #next, plot just the aligned part

    pylab.figure(figsize=(12,12))
    ax = pylab.axes() 
    trans = transforms.blended_transform_factory(ax.transAxes, ax.transData)

    aL = len(aln2fullalnMap)

    #first, get the aligned distance map
    #need map from aln index to structure index
    adistmat = zeros((aL, aL))*nan
    for i,j in [(i,j) for i in range(aL-1) for j in range(i+1,aL)]:
        xi = aln2structMap[i]
        xj = aln2structMap[j]
        if xi == -1 or xj == -1:
            continue
        adistmat[i,j] = distancesM[xi,xj]
        adistmat[j,i] = distancesM[xi,xj]

    #set diagonal for aligned residues
    for i in range(aL):
        if aln2structMap[i] != -1:
            adistmat[i,i] = 0
    acontacts = adistmat < cutoff

    acim = ones(adistmat.shape + (3,))
    acim[acontacts] = ones(3)*0.4
    pylab.imshow(acim, origin='lower', extent=(+0.5,aL+0.5,+0.5,aL+0.5), interpolation='nearest')

    #plot unresolved and deleted regions
    unresolved = where([aln2structMap[n] == -1 for n in range(len(alnseq))])[0]
    fillGroups(unresolved, aL, (0.4,0.4,0.7,0.3))
    deleted = where([c == '-' for c in alnseq])[0]
    fillGroups(deleted, aL, (0.4,0.7,0.4,0.3))
    with open(os.path.join(outpath, 'unresolvedCounts'), 'at') as f:
        print >>f, repr((name, list(unresolved)))
    with open(os.path.join(outpath, 'deleteCounts'), 'at') as f:
        print >>f, repr((name, list(deleted)))

    #plot insertions (lines, red)
    insertpos = searchsorted(aln2fullalnMap, where([x.islower() for x in fullalnseq])[0])
    #for d in set(insertpos):
    #    pylab.axvline(d+0.5, color='r', zorder=9)
    #    pylab.axhline(d+0.5, color='r', zorder=9)
    transY = transforms.blended_transform_factory(ax.transAxes, ax.transData)
    altern = False
    insertlist = [(d,len(list(g))) for d,g in groupby(insertpos)]
    for d, n in insertlist[1:-1]: #first and last are just sequence extension
        pylab.axvline(d+0.5, color='r', zorder=9)
        pylab.axhline(d+0.5, color='r', zorder=9)
        pylab.text(0.005 + altern*0.01, d+0.5, str(n), verticalalignment='center', horizontalalignment='right', transform=transY, path_effects=[PathEffects.withStroke(linewidth=3, foreground="w")], color='r', zorder=100)
        altern = not altern

    with open(os.path.join(outpath, 'insertCounts'), 'at') as f:
        print >>f, repr((name, insertlist))

    #annotate regions
    annotate(regions, lobediv, aL, pylab.gca())
    pylab.title(name + ", aligned sequence only ({}A)".format(cutoff))
    pylab.subplots_adjust(bottom=0.05, right=0.90, top=0.95, left=0.05)
    return adistmat

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

mkdir_p(os.path.join(outpath, 'pdb'))
mkdir_p(os.path.join(outpath, 'aln'))
mkdir_p(os.path.join(outpath, 'adist'))

#erase files
with open(os.path.join(outpath, 'insertCounts'), 'wt') as f:
    pass
with open(os.path.join(outpath, 'deletecounts'), 'wt') as f:
    pass
with open(os.path.join(outpath, 'unresolvedcounts'), 'wt') as f:
    pass

for name in fullalnseqs.keys():
#for name in ['1IEP_A', '2GQG_A']:
    #if os.path.exists(os.path.join(outpath, 'adist/{}.npy'.format(name))):
    #    print name, "Already done"
    #    continue
    if not os.path.exists(os.path.join(distpath, name + '.npy')):
        print name, "Skipping: No coords"
        continue
    print name

    adist = plotContacts(name)
    pylab.subplots_adjust(bottom=0.05, right=0.95, top=0.95, left=0.05)

    pylab.figure(1)
    pylab.savefig(os.path.join(outpath, 'pdb/{}.png'.format(name)))
    pylab.close()
    pylab.figure(2)
    pylab.savefig(os.path.join(outpath, 'aln/{}.png'.format(name)))
    pylab.close()

    pylab.save(os.path.join(outpath, 'adist/{}'.format(name)), adist)
