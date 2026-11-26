"""Dimension removed by the constant-bivector Ward identity."""
from math import comb, prod
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory

import numpy as np

from soft_target import trace_count
from soft_trace_dimensions import extract, multiply_by_slot


def exterior_count(genus, prime):
    n = 2*genus + 4
    shape = tuple(2*genus+int(j < 2) for j in range(5))
    size = tuple(k+5-i for i,k in enumerate(shape))
    assert prime > n and prod(size)*prime*(n+6) < 2**63 and prime**2 < 2**63
    unit = np.zeros(size,dtype=np.int64)
    unit[(0,)*5] = 1
    symmetric = [unit]
    for degree in range(1,n+1):
        total = sum(multiply_by_slot(symmetric[degree-stride],stride,2)
                    for stride in range(1,degree+1))
        symmetric.append((total%prime)*pow(degree,-1,prime)%prime)
    return extract(symmetric[n],shape,prime)


def check():
    with TemporaryDirectory() as tmp:
        exe = Path(tmp)/'character'
        run(['c++','-O3','-std=c++17',str(Path(__file__).with_name('soft_target_char.cpp')),
             '-o',str(exe)],check=True)
        for genus,physical in ((2,3),(3,518),(4,213877),(5,118767557)):
            n = 2*genus+4
            bound = 32*120*prod(comb(n-1+2*genus+4-i+int(i < 2),n-1) for i in range(5))
            assert bound < 2**127
            rank = trace_count(exe,2*genus,n,exterior_rank=2)
            assert rank.denominator==1 and int(rank)%1009 == exterior_count(genus,1009)
            assert int(rank) == {2:0,3:193,4:77199,5:37945455}[genus]
            print(genus, 'constant bivector rank',rank,'reciprocal potentials',physical-rank,flush=True)


if __name__=='__main__':
    check()
