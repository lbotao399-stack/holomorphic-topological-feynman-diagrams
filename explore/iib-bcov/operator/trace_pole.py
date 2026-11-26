"""Residue equations when a pulled-back trace denominator depends on six block sums."""
from itertools import combinations
from random import Random

from flint import fmpz_mat


def check(size):
    # The six nonparallel conserved block sums span exactly four dimensions.
    totals = [tuple(int(i==j) for j in range(5)) for i in range(4)]
    totals.append((1,1,1,1,0))
    totals.append(tuple(-sum(col) for col in zip(*totals)))
    rng = Random(6168+size)
    blocks = []
    for total in totals:
        while True:
            block = [tuple(rng.randint(-5,5) for _ in range(5)) for _ in range(size-1)]
            block.append(tuple(x-sum(col) for x,col in zip(total,zip(*block))))
            if fmpz_mat(block).rank()==min(5,size):
                blocks.append(block)
                break
    edges = list(combinations(range(size),2))
    n = len(blocks)*size
    offset = len(blocks)*len(edges)
    matrix = fmpz_mat(5*n,offset+5+n)
    for b,block in enumerate(blocks):
        for i,row in enumerate(block):
            slot = size*b+i
            for axis in range(5):
                matrix[5*slot+axis,offset+axis] = -1
                matrix[5*slot+axis,offset+5+slot] = -row[axis]
            for e,(j,k) in enumerate(edges):
                if i in (j,k):
                    other = block[k if i==j else j]
                    for axis in range(5):
                        matrix[5*slot+axis,len(edges)*b+e] = other[axis]
    kernel = matrix.ncols()-matrix.rank()
    assert kernel==(6 if size==6 else 0)
    print('Six blocks of',size,'with residue kernel dimension',kernel)


def check_split():
    old = [tuple(int(i==j) for i in range(5)) for j in range(5)]
    summed = tuple(map(sum,zip(*old)))
    q = (1,-1,0,0,0)
    crossings = []
    for t in range(4):
        soft = tuple(t*x for x in q)
        hard = [tuple(x-y for x,y in zip(old[0],soft)),*old[1:]]
        row = [0]*5
        for i,p in enumerate(hard):
            fused = hard[:i]+[tuple(x+y for x,y in zip(p,soft))]+hard[i+1:]
            residue = int(fmpz_mat(fused).det())**2
            if not t:
                assert residue==1
            for j in range(5):
                row[j] += p[j]*residue
        crossings.append(q[0]*row[1]-q[1]*row[0])
    assert crossings[0]==q[0]*summed[1]-q[1]*summed[0]==2
    print('The common fused residue gives soft transverse pairing',crossings[0])


if __name__=='__main__':
    for size in (2,5,6):
        check(size)
    check_split()
