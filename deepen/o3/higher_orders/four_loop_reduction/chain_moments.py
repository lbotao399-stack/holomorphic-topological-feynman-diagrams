"""Integrate all four bivalent-vertex Schwinger ratios of a 5d cubic four-loop graph."""
from fractions import Fraction as F
from functools import cache
from itertools import combinations
from math import factorial, prod
from pathlib import Path
import json
import sys

from flint import fmpq_mpoly_ctx
import sympy as sp

here = Path(__file__).resolve().parent
sys.path.insert(0, str(here.parents[3]/'explore/ht-families/verify'))
import vcore

graphs = json.loads((here.parent/'four-loop-graphs.json').read_text())
variables = fmpq_mpoly_ctx.get(['x0','x1','x2','x3'])
x = variables.gens()


def fraction(q):
    return F(int(q.numerator),int(q.denominator))


def rational(q):
    return variables.constant(int(q.p)) / int(q.q)


def parity(order):
    return (-1)**sum(a>b for i,a in enumerate(order) for b in order[i+1:])


class Reduction:
    def __init__(self, graph):
        self.chains = graph['chains']
        self.core = [(p[0],p[-1]) for p in self.chains]
        edges = [edge for p in self.chains for edge in zip(p,p[1:])]
        self.net = vcore.Net(10,edges)
        self.expanded_cycles = self.net.cycle_basis()
        self.lengths = [len(p)-1 for p in self.chains]
        cursor = 0
        columns = []
        for length in self.lengths:
            column = [row[cursor] for row in self.expanded_cycles]
            assert all([row[i] for row in self.expanded_cycles] == column
                       for i in range(cursor,cursor+length))
            columns.append(column)
            cursor += length
        self.C = sp.Matrix(columns).T
        self.interior = [(e,j,v) for e,p in enumerate(self.chains)
                         for j,v in enumerate(p[1:-1])]
        self.groups = [[i for i,(k,j,v) in enumerate(self.interior) if k==e]
                       for e in range(9)]
        self.core_net = vcore.Net(6,self.core)
        tree,_ = self.core_net.dfs_tree()
        self.tree = sorted(tree)
        self.incidence = sp.Matrix(self.core_net.Bmat()).T
        self.tree_inverse = self.incidence[:,self.tree].inv()
        self.minors = {s:int(self.C[:,s].det()) for s in combinations(range(9),4)}

    def minor(self, cols):
        if len(set(cols)) < 4:
            return 0
        return parity(cols)*self.minors[tuple(sorted(cols))]

    @cache
    def moment(self, powers):
        return F(prod(factorial(n) for n in powers),
                 prod(factorial(len(g)+sum(powers[i] for i in g)) for g in self.groups))

    def integrate_polynomial(self, poly):
        return sum(self.moment(m)*fraction(v) for m,v in poly.to_dict().items())

    def density(self, R, sources, at=None):
        """Return (projective density times sqrt(U), U), before integrating the nine core times."""
        assert len(R)==9 and all(t>0 for t in R)
        assert len(sources)==4 and tuple(map(sum,zip(*sources)))==(0,0)
        C = self.C
        M = C*sp.diag(*R)*C.T
        U = int(M.det()) if all(isinstance(t,int) for t in R) else M.det()
        B = C.T*M.inv()*C
        response = sp.eye(9)-B*sp.diag(*R)
        charges = [[0]*2 for _ in self.core]
        offsets = [[variables.constant(0) for _ in range(2)] for e in self.core]
        vertical = []
        for i,(e,j,v) in enumerate(self.interior):
            p = sources[v-6]
            charges[e] = [a+b for a,b in zip(charges[e],p)]
            vertical.append((e,tuple(charges[e])))
            for k in range(2):
                offsets[e][k] += x[i]*charges[e][k]
        for k in range(2):
            injections = [sum(charges[e][k] for e,(_,v) in enumerate(self.core) if v==j)
                          for j in range(5)]
            for e,n in zip(self.tree,self.tree_inverse*sp.Matrix(injections)):
                offsets[e][k] += int(n)
        horizontal = [[sum((rational(response[e,j])*offsets[j][k] for j in range(9)),
                            variables.constant(0)) for k in range(2)] for e in range(9)]
        # First eight horizontal R directions, then the four vertical x directions.
        result = variables.constant(0)
        for J in combinations(range(8),4):
            d = self.minors[J]
            if not d:
                continue
            a,b,c,f = J
            haf = B[a,b]*B[c,f]+B[a,c]*B[b,f]+B[a,f]*B[b,c]
            if not haf:
                continue
            H = [j for j in range(8) if j not in J]
            labels = H + [e for e,_ in vertical]
            coefficient = [0]*16
            for S in combinations(range(8),4):
                T = tuple(i for i in range(8) if i not in S)
                minor = self.minor([labels[i] for i in S])*self.minor([labels[i] for i in T])
                if minor:
                    bits = sum(1<<i for i in T if i<4)
                    coefficient[bits] += parity(S+T)*minor*prod(
                        momentum[int(4+j in T)] for j,(e,momentum) in enumerate(vertical))
            poly = sum((a*prod(horizontal[e][(m>>j)&1] for j,e in enumerate(H))
                        for m,a in enumerate(coefficient) if a), variables.constant(0))
            result += rational(d*haf*parity(list(J)+H)/4)*poly
        # Source chart interleaves each chain's interior splits with its chain time.
        # Both the chart orientation and the absolute pushforward are retained by this sign.
        sign = (-1)**sum((8-e)*len(g) for e,g in enumerate(self.groups))
        factor = F(sign*prod(R[e] for e,_ in vertical),1)/(F(R[-1])*F(U)**2)
        val = self.integrate_polynomial(result) if at is None else fraction(result(*at))
        return factor*val, U


def verify():
    from flint import fmpq
    sources = [(2,1),(1,4),(-2,3),(-1,-8)]
    R = [2,3,5,7,11,13,17,19,23]
    for idx in (0,1,2,3,11,39,49):
        red = Reduction(graphs[idx])
        at = [F(i+1,17) for i in range(4)]
        times = []
        for e,g in enumerate(red.groups):
            shares = [at[i] for i in g]
            times.extend(R[e]*v for v in [1-sum(shares)]+shares)
        expected,U,_ = vcore.density_slice(red.net,red.expanded_cycles,times,
                                           [(0,0)]*6+sources[:3])
        jacobian = prod(R[e]**len(group) for e,group in enumerate(red.groups))
        actual,base_U = red.density(R,sources,at=[fmpq(a.numerator,a.denominator) for a in at])
        assert U == base_U and actual == expected*jacobian, (idx,actual,expected*jacobian)
        val,_ = red.density(R,sources)
        print(idx, val, 'sqrt(U) times eight-dimensional density, U =',U)


if __name__=='__main__':
    verify()
