#!/usr/bin/env python
from pymol.cgo import *
from pymol import cmd
import sys,re
from scipy import array, zeros
import numpy as np
from pdbParse import loadPDB

# pymol scripy to show cylinders for interactions.

colors = {}
colors['red']   = [1.0,0.0,0.0]
colors['green'] = [0.0,1.0,0.0]
colors['blue']  = [0.0,0.0,1.0]
colors['cyan']  = [0.5,0.5,1.0]

interactions = np.load('wFB_kinase.npy')
L = int(((1+np.sqrt(1+8*len(interactions)))//2) + 0.5)

with open('../structure/03-21-2019/alnMap', 'rb') as f:
    lines = [l.strip(b'\n').split(b'|') for l in f]
    m = dict([(l[0].decode('ascii'), l[1:]) for l in lines])

pdbname = sys.argv[1]
chain = sys.argv[2]

chainsel = 'chain {}'.format(chain)

cmd.load(pdbname)
pdb = loadPDB(pdbname)
pdb = pdb[pdb.chain == chain.encode('ascii')]

CB = ((pdb.atom == b' CB ') |
      ((pdb.atom == b' CA ') & (pdb.resName == b'GLY')))

structname = os.path.splitext(os.path.basename(pdbname))[0]
amap = m[structname+'_'+chain]

cColor = colors['cyan']

# show potts interactions as sticks
cyls = [COLOR,0.0,0.0,0.0,]
pairs = ((i,j) for i in range(L-1) for j in range(i+1,L))
for n,(i,j) in enumerate(pairs):
    p = 1.5*(interactions[n] - 0.14)
    if p < 0:
        continue
    a1 = pdb[(pdb.resID == amap[i]) & CB]
    a2 = pdb[(pdb.resID == amap[j]) & CB]
    if len(a1) < 1 or len(a2) < 1:
        continue
    cyls.extend([CYLINDER] +  list(a1[0].coords) + list(a2[0].coords) + 
                [p] + cColor + cColor)
cmd.load_cgo(cyls, 'Potts')

# show conservation wit spheres
qeff = np.load('qeff.npy')
spheres = [COLOR,0.0,0.0,0.0,]
for i,q in enumerate(qeff):
    a = pdb[(pdb.resID == amap[i]) & CB]
    if len(a) < 1:
        continue
    spheres.extend([SPHERE] + list(a[0].coords) + [1.0/q])
cmd.load_cgo(spheres, 'Cnsrv')

# label residue numbers
for i in range(L):
    mi = amap[i].decode('ascii')
    if mi.strip() != '-':
        sele = 'resi {} and (name CB or (name CA and resn GLY))'.format(mi)
        sele = chainsel + " and " + sele
        cmd.label(sele, '"{}"'.format(str(i)))


# color regions of kinase
cmd.bg_color('white')
cmd.hide('cartoon')
cmd.hide('sticks', 'not {}'.format(chainsel))
cmd.hide('spheres', 'not {}'.format(chainsel))

regions = [
        ('N Lobe',     (  0,  67), 'palegreen'),
        ('C Lobe',     ( 67, 175), 'white'),
        ('P-loop',     (  0,   9), 'yelloworange'),
        ('b-3 K',      ( 23,  24), 'skyblue'),
        ('alpha-C',    ( 30,  45), 'orange'),
        ('alpha-C E',  ( 37,  38), 'deeppurple'),
        ('Gatekeeper', ( 66,  67), 'black'),
        ('Hinge',      ( 67,  74), 'gray50'),
        ('HRD',        (108, 111), 'magenta'),
        ('cata',       (111, 117), 'yellow'),
        ('cata N',     (115, 116), 'chartreuse'),
        ('cata2',      (117, 128), 'cyan'),
        ('DFG',        (128, 131), 'red'),
        ('Act Loop',   (131, 151), 'blue'),
        ('APE',        (151, 154), 'deepteal')
    ]

for name, (s, e), color in regions:
    rs = '-'
    while rs == '-' and s < e:
        rs = amap[s].decode('ascii').strip()
        s += 1

    re = '-'
    while re == '-' and e > s:
        re = amap[e-1].decode('ascii').strip()
        e -= 1

    print(name, s, e)
    if re == '-':
        continue

    sele = 'resi {}-{}'.format(rs, re)
    sele = chainsel + " and " + sele
    cmd.show('cartoon', sele)
    cmd.color(color, sele)

cmd.center(chainsel)
cmd.hide('(solvent)')
