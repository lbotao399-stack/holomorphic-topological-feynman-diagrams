"""Characters of quadratic-soft genus-three bulk potentials and their gauge shifts."""
from collections import Counter
from fractions import Fraction
from math import comb,factorial,prod
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory
from genus_two_probes import partitions


def multiplicities(exe,weight=6,n=10,cut=2):
    query=[]
    weights=[]
    for cycle in partitions(n):
        order=prod(m**f*factorial(f) for m,f in Counter(cycle).items())
        query.append((weight,1,1,len(cycle),*(z for j in cycle for z in (j,cut))))
        weights.append((Fraction(cycle.count(1),order),Fraction(1,order),Fraction()))
    for cycle in partitions(n-1):
        order=prod(m**f*factorial(f) for m,f in Counter(cycle).items())
        query.append((weight,0,1,len(cycle)+1,1,max(0,cut-1),
                      *(z for j in cycle for z in (j,cut))))
        weights.append((Fraction(),Fraction(),Fraction(1,order)))
    res=run([exe],text=True,input='\n'.join(' '.join(map(str,q)) for q in query),
            capture_output=True,check=True)
    numbers=list(map(int,res.stdout.split()))
    return tuple(sum(m*c[i] for m,c in zip(numbers,weights)) for i in range(3))


def trace_count(exe,weight,n,exterior_rank=0):
    query=[]
    weights=[]
    for cycle in partitions(n):
        order=prod(m**f*factorial(f) for m,f in Counter(cycle).items())
        query.append((weight,exterior_rank,1,len(cycle),*(z for length in cycle for z in (length,2))))
        weights.append(Fraction(1,order))
    result=run([exe],text=True,input='\n'.join(' '.join(map(str,row)) for row in query),
               capture_output=True,check=True)
    return sum(v*w for v,w in zip(map(int,result.stdout.split()),weights))


def genus_counts(exe,g):
    # The complete linear coefficient is bounded by the uncensored tensor
    # character. Five momentum relations and the Weyl antisymmetrizer add
    # factors at most 32 and 120 to this absolute bound.
    n=2*g+4
    bound=32*120*prod(comb(n-1+2*g+4-i+int(i==0),n-1) for i in range(5))
    if bound>=2**127:
        raise OverflowError("Character exceeds the proven signed-128-bit bound.")
    trace=trace_count(exe,2*g,n-1)
    array,common,longitudinal=multiplicities(exe,2*g,n,2)
    return trace,array,common,longitudinal,array-common-longitudinal


if __name__=='__main__':
    with TemporaryDirectory() as tmp:
        exe=Path(tmp)/'soft-target-character'
        run(['c++','-O3','-std=c++17',str(Path(__file__).with_name('soft_target_char.cpp')),
             '-o',str(exe)],check=True)
        for cut in (0,1,2):
            vals=multiplicities(exe,cut=cut)
            assert vals==((4013,493,795),(2774,335,532),(770,96,156))[cut]
            print('Minimum momentum degree',cut,'potential, common, longitudinal',vals,'quotient',
                  vals[0]-sum(vals[1:]),flush=True)
