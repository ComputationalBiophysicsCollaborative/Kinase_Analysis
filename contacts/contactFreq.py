#!/usr/bin/env python
# coding=utf-8
from scipy import *
import scipy.io as sio
import matplotlib
import pylab
import matplotlib.patheffects as PathEffects
import matplotlib.transforms as transforms
import numpy as np
import sys, os, re, errno
from annotate import annotate
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seqtools
from Bio.Alphabet import IUPAC
alpha = '-' + IUPAC.protein.letters

matplotlib.rc('savefig', dpi=200)

pdblistfn = sys.argv[1]
pdbseqfn = sys.argv[2]
weight = float(sys.argv[3])
regionfn = sys.argv[4]
distpath = sys.argv[5]
outpath = sys.argv[6]
L = int(sys.argv[7])

def finishfigs():
    #pylab.show()
    pylab.close('all')

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

mkdir_p(outpath)

with open(pdblistfn) as f:
    pdblist = [l.strip() for l in f]

with open(pdbseqfn) as f:
    pdbseq = dict(l.split() for l in f)

with open(regionfn) as f:
    dat = [l.split(',') for l in f.readlines()]
    regions = [(n,int(s),int(e)) for n,s,e in dat]
lobediv = [e for n,s,e in regions if n == 'Hinge'][0]

#with open(os.path.join(distpath, 'insertCounts')) as f:
#    insertdat = [eval(l) for l in f.readlines()]
#    insertdat = {name: array(d) for name, d in insertdat}
#    inserts = zeros((len(insertdat), L), dtype=int)
#    for n,(name,d) in enumerate(insertdat):
#        d = array(d)
#        inserts[n,d[:,0]-1] = d[:,1]

#with open(os.path.join(distpath, 'deleteCounts')) as f:
#    deletedat = [eval(l) for l in f.readlines()]
#    deletes = zeros((len(deletedat), L), dtype=int)
#    for n,(name,d) in enumerate(deletedat):
#        deletes[n,array(d, dtype=int)-1] += 1

#with open(os.path.join(distpath, 'unresolvedCounts')) as f:
#    unresolveddat = [eval(l) for l in f.readlines()]
#    unresolveds = zeros((len(unresolveddat), L), dtype=int)
#    for n,(name,d) in enumerate(unresolveddat):
#        unresolveds[n,array(d, dtype=int)-1] += 1
##assume all the count files have same order
#freqnames = [n for n,d in unresolveddat] 

distances = []
structs = []
for struct in pdblist:
    try:
        dist = load(os.path.join(distpath, "adist", struct+".npy"))
    except:
        print("no distances for", struct)
        continue
    distances.append(dist)
    structs.append(struct)
seqs = [pdbseq[s] for s in structs]
seqs = np.array([[alpha.index(c) for c in s] for s in seqs], dtype='u1')
weights = 1.0/seqtools.nsim(seqs, int(L*(weight)))
print(weight, sum(weights))
distances = array(distances).transpose((1,2,0))
valid = ~isnan(distances)

ntn = nan_to_num 
nrm = dot(valid, weights)
totmean = dot(ntn(distances), weights)/nrm
totstd = sqrt(dot(ntn(distances - totmean[:,:,None])**2, weights)/nrm)

def plotMap(im, vmin=0,vmax=1, cmap='bwr', nocbar=False, title=None,
            freqmask=slice(0,2**63-1)):
    #see http://matplotlib.org/1.4.2/mpl_toolkits/axes_grid/users/overview.html
    pylab.figure(figsize=(16,12))

    axmap = pylab.subplot(111)
    divider = make_axes_locatable(axmap)
    if not nocbar:
        cax = divider.append_axes("right", size="5%", pad=0.5)
    axins = divider.append_axes("left", size=0.5, pad=0.5, sharey=axmap)
    axdel = divider.append_axes("left", size=0.5, pad=0.15, sharey=axmap)
    axunr = divider.append_axes("left", size=0.5, pad=0.15, sharey=axmap)

    im = axmap.imshow(im, origin='lower', extent=(+0.5,L+0.5,+0.5,L+0.5), 
                      interpolation='nearest', cmap=cmap, 
                      vmin=vmin, vmax=vmax, aspect='equal')
    annotate(regions, lobediv, L, axmap)
    axmap.set_aspect(1)
    
    cbar = None
    if not nocbar:
        cbar = pylab.colorbar(im, cax=cax)
    if title:
        axmap.set_title(title)
    
    def wmean(s):
        num = sum((s*weights[:,newaxis])[freqmask], axis=0)
        return num/sum(freqweights[freqmask])

    #axins.barh(arange(L)+1, wmean(inserts), height=1, color='gray', ec='none')
    #axins.set_xlim(0,6)
    #axins.set_xticks([0,3,6])
    #axins.invert_xaxis()
    #axins.text(6*0.9, 1, 'Mean Number of Insertions', rotation='vertical',
    #           verticalalignment='bottom')
    #pylab.setp(axins.get_yticklabels(), visible=False)


    #axdel.barh(arange(L)+1, wmean(deletes), height=1, color='gray', ec='none')
    #axdel.set_xlim(0,1)
    #axdel.set_xticks([0,1])
    #axdel.invert_xaxis()
    #axdel.text(0.9, 1, 'Frequency of Deletions', rotation='vertical',
    #           verticalalignment='bottom')
    #pylab.setp(axdel.get_yticklabels(), visible=False)

    #axunr.barh(arange(L)+1, wmean(unresolveds), height=1,
    #           color='gray', ec='none')
    #axunr.set_xlim(0,1)
    #axunr.set_xticks([0,1])
    #axunr.invert_xaxis()
    #axunr.text(0.9, 1, 'Frequency of Unresolved', rotation='vertical',
    #           verticalalignment='bottom')
    #pylab.setp(axunr.get_yticklabels(), visible=False)

    #axmap.set_yticks(axmap.get_xticks())
    axmap.set_xlim(0.5,L+0.5)
    axmap.set_ylim(0.5,L+0.5)
    return axmap, axins, axdel, cbar

pylab.figure()

distv = totmean[~isnan(totmean)]
plotMap(totmean, vmin=min(distv), vmax=max(distv), cmap='YlGnBu_r', title='Mean Distance (All Seqs)')
pylab.savefig(os.path.join(outpath, 'tot_dist.png'))
save(os.path.join(outpath, 'totdist'), totmean)

distv = totstd[~isnan(totmean)]
plotMap(totstd, vmin=min(distv), vmax=max(distv), cmap='YlGnBu',
        title='Variance in Distance (All Seqs)')
pylab.savefig(os.path.join(outpath, 'tot_var.png'))

distv = (totstd/totmean)[~isnan(totmean)]
plotMap(totstd/totmean, vmin=min(distv), vmax=max(distv), cmap='YlGnBu',
        title='Variance/Mean in Distance (All Seqs)')
pylab.savefig(os.path.join(outpath, 'tot_var_normed.png'))

finishfigs()

for cutoff in [5.0, 5.25, 5.5, 6.0, 6.5, 7.0, 8.0, 10.0]:
    totp = dot(distances < cutoff, weights)/nrm
    
    cutoffdir = os.path.join(outpath, str(cutoff)+'A')
    mkdir_p(cutoffdir)

    plotMap(totp, cmap='gray_r', 
            title='Total Contact Frequency, {} Å cutoff'.format(cutoff))
    pylab.savefig(os.path.join(cutoffdir, 'contacts_tot.png'))
    save(os.path.join(cutoffdir, 'totcontactfreq'), totp)
    pylab.close()


