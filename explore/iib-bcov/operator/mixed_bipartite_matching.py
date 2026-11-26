"""Integrate the ordered positions of the mixed colors on a bipartite perfect matching."""
from argparse import ArgumentParser
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import islice, permutations
from math import factorial, prod
from pathlib import Path
from subprocess import run
import gzip
import os
import gc
import sys

import numpy as np
from flint import fmpq
from mixed_bare_symbolic import classified_orbits
from mixed_bipartite_one_star import METRIC, PRIME, PrimitiveChamber, residue
from mixed_prism_one_star import chamber_rows
from mixed_prism_rung_sum import ColorMoments


def matrix_orderings(core):
    sides=((0,4,5),(1,2,3))
    return [[core.index(tuple(sorted((i,j)))) for i in row for j in col]
            for left,right in (sides,sides[::-1])
            for row in permutations(left) for col in permutations(right)]


def project(marks,destination,color_word=None):
    probe=PrimitiveChamber(4,marks)
    orderings=matrix_orderings(probe.core)
    powers16=np.array([16**i for i in reversed(range(9))],dtype=np.int64)
    adapted=os.environ.get('MIXED_ADAPTED_CYCLES','1' if len(set(marks))==3 else '0')=='1'
    batch=ColorMoments(marks,(),4,adapted=adapted)
    checks=[(m,1,1,tuple(map(fmpq,m)),[0,Fraction()]) for m in (METRIC,(11,5,2,7,9,6,8,3,4))]
    polynomials={}
    def collect(part):
        for degree,poly in zip(part.coefficients,part.projected(),strict=True):
            if degree in polynomials:polynomials[degree]+=poly
            else:polynomials[degree]=poly
    rows=list(chamber_rows(marks,4))
    if color_word is not None and color_word not in {color for color,_ in rows}:
        raise ValueError("color is not a canonical chamber of this placement")
    for color,count in rows:
        if color_word is not None and color!=color_word:continue
        batch.add(color,Fraction(count,len(probe.actions)),checks,collect)
    expected=[r for r in classified_orbits(4)[0] if r[0]==marks]
    wanted={tuple(r[1]) for r in expected}
    if color_word is not None:wanted.intersection_update((color_word,))
    assert wanted==set(batch.nonzero)
    print('Hard colors passing the independent census:',len(wanted),flush=True)
    weights=defaultdict(fmpq)
    count=0
    for poly in polynomials.values():
        for _,_,_,metric,state in checks:
            state[1]+=Fraction(str(poly(*metric)))
        stream=iter(poly.terms())
        while group:=list(islice(stream,8192)):
            ex=np.array([list(m) for m,_ in group],dtype=np.int64)
            assert (ex.sum(axis=1)==43).all() and ex.min()>=0 and ex.max()<12
            powers=ex+1
            keys=np.minimum.reduce([powers[:,o]@powers16 for o in orderings])
            for key,(_,coefficient) in zip(keys.tolist(),group,strict=True):
                weights[key]+=coefficient
            count+=len(group)
    for _,_,_,_,state in checks:assert state[0]==residue(state[1])
    stats=Counter()
    with gzip.open(destination.with_suffix('.partial'),'wt') as output:
        for key,coefficient in sorted(weights.items()):
            a=tuple((key>>(4*i))&15 for i in reversed(range(9)))
            rows=[sum(a[i:i+3]) for i in (0,3,6)]
            cols=[sum(a[i::3]) for i in range(3)]
            assert max(rows+cols)<26 and all(52-rows[i]-cols[j]+a[3*i+j]>13
                    for i in range(3) for j in range(3))
            if not coefficient:continue
            stats[min(rows+cols)]+=1
            norm=fmpq(prod(factorial(i-1) for i in a),factorial(12))
            print(key,norm*coefficient,file=output)
    destination.with_suffix('.partial').replace(destination)
    print('Raw monomials and nonzero core orbits by smallest vertex sum:',count,sorted(stats.items()),flush=True)


def collect_colors(word,destination):
    marks=tuple(map(int,word))
    weights=defaultdict(fmpq)
    selected=[r for r in classified_orbits(4)[0] if r[0]==marks]
    for _,color,_,_ in selected:
        label=''.join(map(str,color))
        file=destination.with_name(f'{word}-{label}.gz')
        if not file.exists():
            run([sys.executable,'-B','-u',__file__,word,str(file),'--color',label],check=True)
        with gzip.open(file,'rt') as source:
            for row in source:
                key,weight=row.split()
                weights[int(key)]+=fmpq(weight)
    pending=destination.with_suffix('.partial')
    with gzip.open(pending,'wt') as out:
        for k,w in sorted(weights.items()):
            if w:print(k,w,file=out)
    pending.replace(destination)
    print('Nonzero colors and combined core orbits:',len(selected),sum(bool(v) for v in weights.values()),flush=True)


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('word',help='six sorted core path labels of a canonical occupancy')
    parser.add_argument('destination',type=Path)
    parser.add_argument('--color',help='project just one hard-color orbit')
    args=parser.parse_args()
    if len(args.word)!=6 or any(c not in '012345678' for c in args.word) or list(args.word)!=sorted(args.word):
        parser.error('a placement is a sorted six-digit path word')
    if tuple(map(int,args.word)) not in {r[0] for r in classified_orbits(4)[0]}:
        parser.error('the placement has no surviving canonical chamber')
    if args.color:
        # Every color is computed in an isolated process. The moment tables are
        # acyclic; scanning their millions of tuples during contraction dominates
        # runtime, while reference counting still releases each collected root.
        gc.disable()
        project(tuple(map(int,args.word)),args.destination,tuple(map(int,args.color)))
    else:
        collect_colors(args.word,args.destination)
