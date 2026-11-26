"""Certify all-momentum tests for one vector and eight scalar slots."""
from fractions import Fraction
from math import factorial, prod
from collections import Counter
from pathlib import Path
from random import Random
from subprocess import run
from tempfile import TemporaryDirectory

from genus_two_probes import det, partitions


def nine(a, b, c):
    p = [tuple(int(i==j) for i in range(5)) for j in range(5)] + [a,b,c]
    return tuple(p + [tuple(-sum(z) for z in zip(*p))])


FIXTURES = (
    nine((1,-1,1,0,-1),(0,1,1,-1,-1),(0,1,-1,2,1)),
    nine((1,2,-1,0,-1),(-2,0,1,-1,2),(0,1,2,-2,1)),
    nine((2,-1,0,-2,1),(0,2,1,-1,-2),(1,-2,2,1,1)))


def pivot_rows(rows):
    basis = {}
    ids = []
    for i,r in enumerate(rows):
        r = [x%1009 for x in r]
        for j in sorted(basis):
            a = r[j]
            r = [(x-a*y)%1009 for x,y in zip(r,basis[j])]
        if any(r):
            j = next(i for i,x in enumerate(r) if x)
            a = pow(r[j],-1,1009)
            basis[j] = [x*a%1009 for x in r]
            ids.append(i)
    return ids


def multiplicities(n,shape):
    t = a = Fraction()
    for cycle in partitions(n):
        h = [1,-1]+[0]*(max(shape)+len(shape)-2)
        for m in cycle:
            for j in range(m,len(h)):
                h[j] += h[j-m]
        ch = det([[h[v-i+j] if v-i+j>=0 else 0 for j in range(5)] for i,v in enumerate(shape)])
        order = prod(m**f*factorial(f) for m,f in Counter(cycle).items())
        t += Fraction(ch,order)
        a += Fraction(ch*cycle.count(1),order)
    return t,a


def columns(executable,p,orbits):
    values = [*p,(len(orbits),),*orbits]
    stdin = '\n'.join(' '.join(map(str,row)) for row in values)
    result = run([executable],input=stdin,text=True,capture_output=True,check=True)
    return [tuple(map(int,line.split())) for line in result.stdout.splitlines()]


def check():
    assert multiplicities(9,(5,4,4,4,4)) == (13,75)
    assert multiplicities(9,(4,4,4,4,4)) == (7,24)
    assert multiplicities(8,(4,4,4,4,4)) == (4,12)
    rng = Random(93576)
    orbits = [tuple(sorted(rng.randrange(126) for j in range(4))) for i in range(640)]
    for _ in range(96):
        a,b,c = [rng.randrange(126) for j in range(3)]
        orbits += [(a,a,a,a),(a,a,b,b),(a,a,b,c),(a,a,a,b)]
    selected_columns = (1,2,3,4,6,7,8,9,10,11,13,15,18,23,25,27,28,29,30,32,33,34,35,36,37,39,41,42,47,48,51,53,69,78,82,85,86,91,105,115,148,158,159,178,232,353,371,427,479,583,640)
    orbits = [orbits[i] for i in selected_columns]
    with TemporaryDirectory() as tmp:
        exe = Path(tmp)/'quartic-probe-orbits'
        source = Path(__file__).with_name('quartic_probe_orbits.cpp')
        run(['c++','-O3','-std=c++17',str(source),'-o',str(exe)],check=True)
        rows = []
        for p in FIXTURES:
            rows.extend(zip(*columns(exe,p,orbits)))
            ids = pivot_rows(rows)
            print('Rank after',len(rows)//36,'fixtures:',len(ids),flush=True)
            if len(ids)==51:
                break
    assert len(ids)==51
    selected = [rows[i] for i in ids]
    cols = pivot_rows(zip(*selected))
    residue = det([[r[j] for j in cols] for r in selected])%1009
    assert residue==714
    print('The exact minor is 714 modulo 1009.')


if __name__=='__main__':
    check()
