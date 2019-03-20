#!/usr/bin/env python2
from scipy import *
import sys, re
import pylab

def matchstr(s):
    m = list(re.finditer('[A-Z-]', s))
    return s[m[0].start():m[-1].end()]

rawseqs = []
with open(sys.argv[1]) as f:
    while True:
        name = f.readline()
        seq = f.readline()
        if name == '':
            break
        rawseqs.append((name.split()[0][1:], seq.strip()))

dat = []
alignedseqs = []
fullseqs = []
for name, seq in rawseqs:
    noends = matchstr(seq)
    aligned = re.sub('[a-z]', '', noends)
    #print "".join(aligned[i] for i in [113,114,115,133,134,135]), aligned[133]
    if re.search('[XZB]', aligned):
        continue
    if any([aligned[i] == '-' for i in [113,114,115,133,134,135]]):
        continue
    if aligned[133] != 'D':
        continue
    if aligned.count('-') > 10:
        continue
    if len(noends) - len(aligned) > 40:
        continue
    #print ">{}\n{}".format(name, noends)
    #print ">{}\n{}".format(name, aligned)
    alignedseqs.append((name, aligned))
    fullseqs.append((name, noends))
    dat.append((aligned.count('-'), len(noends) - len(aligned)))
dat = array(dat)
#print dat.shape

with open(sys.argv[2], 'wt') as f:
    f.write("\n".join([">{}\n{}".format(*s) for s in alignedseqs]))

with open(sys.argv[3], 'wt') as f:
    f.write("\n".join(["{} {}".format(*s) for s in fullseqs]))

with open(sys.argv[4], 'wt') as f:
    f.write("\n".join([s[1] for s in alignedseqs]))

fig = pylab.figure(figsize=(16,8))
pylab.subplot(221)
pylab.hist(dat[:,0], bins=arange(0, max(dat[:,0])))
pylab.xlabel('Deletions')
pylab.subplot(223)
pylab.hist(dat[:,1], bins=arange(0, max(dat[:,1])))
pylab.xlabel('Insertions')
pylab.subplot(122)
p = max(dat.flatten())
pylab.hist2d(dat[:,0], dat[:,1], bins=arange(p+1), cmap='gray_r')
pylab.xlabel('Deletions')
pylab.ylabel('Insertions')
pylab.savefig('ins-dels_dist.png')
pylab.show()
