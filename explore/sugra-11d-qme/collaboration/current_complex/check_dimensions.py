"""Compare finite-torus characters to direct jet-monomial counting."""
from collections import defaultdict
from itertools import permutations
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'cohomology'))
from jets import multi_indices
from relative_dimensions import dimensions

def direct(n,degree,qmax,k):
    zero=(0,)*n
    counts=[[defaultdict(int) for _ in range(degree+1)] for _ in range(qmax+1)]
    counts[0][0][zero]=1
    for h in range(1,degree+1):
        for alpha in multi_indices(n,h+1):
            for axis in range(n):
                if axis==0 and alpha[0]:
                    continue
                weight=tuple(int(i==axis)-a for i,a in enumerate(alpha))
                for q in range(qmax,0,-1):
                    for d in range(degree,h-1,-1):
                        dest=counts[q][d]
                        for t,v in counts[q-1][d-h].items():
                            dest[tuple(a+b for a,b in zip(t,weight))]+=v
    ans=[[0]*(qmax+1) for _ in range(2)]
    for p in permutations(range(n)):
        sign=(-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
        t=tuple(k+p[i]-i for i in range(n))
        for q in range(qmax+1):
            ans[1][q]+=sign*counts[q][degree].get(t,0)
            ans[0][q]+=sign*sum(counts[q][degree-1].get(tuple(a+int(i==j) for i,a in enumerate(t)),0)
                                 for j in range(n))
    return ans

if __name__=='__main__':
    for n,k,qmax in [(2,-3,5),(3,-2,5),(5,-1,4)]:
        a=dimensions(n,32,qmax,k)
        assert a==direct(n,-k*n,qmax,k)
        print(n,k,a)
