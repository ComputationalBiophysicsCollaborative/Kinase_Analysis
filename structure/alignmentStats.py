#!/usr/bin/env python2
from scipy import *
import pylab
import re, sys

def matchstr(s):
    m = list(re.finditer('[A-Z-]', s))
    return s[m[0].start():m[-1].end()]

with open(sys.argv[1]) as f:
    seqs = [l.split() for l in f.readlines()]

if len(sys.argv) > 2:
    weights = load(sys.argv[2])
else:
    weights = None

sub = [matchstr(v[1]) for v in seqs]

deletes = array([len(re.findall('-', s)) for s in sub])
inserts = array([len(re.findall('[a-z]', s)) for s in sub])
print inserts

fig = pylab.figure()
fig.suptitle('Weighted Counts')
pylab.subplot(211)
pylab.hist(inserts, bins=30, weights=weights)
pylab.xlabel('Insertions')
pylab.subplot(212)
pylab.hist(deletes, bins=30, weights=weights)
pylab.xlabel('Deletions')

fig = pylab.figure()
fig.suptitle('Raw Counts')
pylab.subplot(211)
pylab.hist(inserts, bins=30)
pylab.xlabel('Insertions')
pylab.subplot(212)
pylab.hist(deletes, bins=30)
pylab.xlabel('Deletions')

pylab.figure()
pylab.plot(inserts, deletes, '.')
pylab.xlabel('Insertions')
pylab.ylabel('Deletions')

pylab.show()
