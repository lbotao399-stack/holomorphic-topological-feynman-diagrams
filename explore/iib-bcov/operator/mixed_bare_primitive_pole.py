"""Kirchhoff pole of two position-integrated primitive mixed Ward chambers."""

from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, permutations, product
from math import factorial, prod
from pathlib import Path
import argparse
import json

from flint import nmod_mat, nmod_mpoly_ctx
from fourfold_quartic_orbits import oriented_symmetries
from fourfold_mixed_soft_ward import V, integrated_soft_placements
from mixed_bare_ranks import cubic_cores, hard_waves, WAVES
from quartic_seven_plane_double import G, fixture, forest_data, image, numerator, solve_graph, tree_flow


PRIME = 1009
COLOR = (0,0,2,1,1,2)
METRIC = (2,3,5,7,11,13,17,19,23)


def residue(value):
    return value.numerator*pow(value.denominator,-1,PRIME)%PRIME


@lru_cache(maxsize=3)
def hard_resistances(edges):
    return forest_data(edges,12)[2]


class PrimitiveChamber:
    def __init__(self,index,marks=None,color=COLOR):
        self.color = color
        self.core = cubic_cores()[index]
        self.marks = marks if marks is not None else ((0,0,0,6,6,6) if index == 3 else (0,0,0,5,5,5))
        self.cycles = G(6,self.core).C
        self.edges = fixture(self.core,self.marks)
        self.graph = G(12,self.edges)
        self.graph.C = [[r[e] for e in range(9) for _ in range(1+self.marks.count(e))] for r in self.cycles]
        (u,t,z),waves = hard_waves(self.color)
        self.source = [(0,)*4]*6 + waves.tolist()
        self.sheared = [(a,b+c,d,c) for a,b,c,d in self.source]
        self.flow = tree_flow(self.graph,self.sheared)
        resistance = hard_resistances(tuple(self.edges))
        def cross(x,y):
            return sum(resistance[6+i][6+j] for i,j in product(x,y))
        four_voltages = 4*resistance[6+u[0]][6+u[1]]+cross(t,z)-cross(u,z)-cross(t,u)
        self.voltage_terms = [(tuple(i for i,p in enumerate(ex) if p),int(k)) for ex,k in four_voltages.terms()]
        assert all(len(ex)==5 for ex,_ in self.voltage_terms)
        self.actions = tuple(oriented_symmetries(self.core))
        self.ctx = nmod_mpoly_ctx.get(tuple(f'x{i}' for i in range(6)),PRIME)
        self.groups = [tuple(i for i,p in enumerate(self.marks) if p==e) for e in range(9)]
        config = tuple(tuple(self.color[i] for i in group) for group in self.groups)
        self.weight = F(8*len({image(config,g) for g in self.actions}),len(self.actions))


    def kirchhoff(self,times):
        return int(nmod_mat([[sum(a*b*t for a,b,t in zip(i,j,times)) for i in self.cycles]
                             for j in self.cycles],PRIME).det())

    def divisor(self):
        times = list(METRIC)
        times[0] = 0
        intercept = self.kirchhoff(times)
        times[0] = 1
        slope = self.kirchhoff(times)-intercept
        times[0] = -intercept*pow(slope,-1,PRIME)%PRIME
        assert not self.kirchhoff(times)
        return tuple(times)

    def direct(self,times,square,voltage):
        points = tuple(F(j+1,len(group)+1) for group in self.groups for j in range(len(group)))
        lengths = []
        for e,group in enumerate(self.groups):
            ends = (F(),*[points[i] for i in group],F(1))
            lengths.extend(times[e]*(b-a) for a,b in zip(ends,ends[1:]))
        flow,_ = solve_graph(self.edges,self.source,lengths,0)
        electrical = integrated_soft_placements(flow,self.graph.C,lengths,V)
        tree = self.kirchhoff(times)
        at = [residue(p) for p in points]
        assert 2*int(square(*at))%PRIME == residue(electrical)*pow(tree,8,PRIME)%PRIME
        (u,t,z),_ = hard_waves(self.color)
        signed = F()
        for (a,b),(c,d),(f,g) in product(permutations(u),permutations(t),permutations(z)):
            source = self.source.copy()
            for site,wave in zip((a,c,d,f,g,b),WAVES[:6]):
                source[6+site] = wave
            _,phi = solve_graph(self.edges,source,lengths,0)
            signed += phi[6+a][3]-phi[6+c][3]
        assert 2*int(voltage(*at))%PRIME == residue(signed)*tree%PRIME

    def numerator(self,times,check=False):
        x = self.ctx.gens()
        segments = []
        for e,group in enumerate(self.groups):
            ends = (self.ctx.constant(0),*[x[i] for i in group],self.ctx.constant(1))
            segments.extend(times[e]*(b-a) for a,b in zip(ends,ends[1:]))
        minors = {(i,j):numerator(self.graph,self.sheared,self.flow,[(3,i),(3,j)],edge_times=segments)
                  for i,j in combinations(range(4),2)}
        square = self.ctx.constant(0)
        for p,q in combinations(range(9),2):
            minor = sum((self.cycles[i][p]*self.cycles[j][q]-self.cycles[i][q]*self.cycles[j][p])*d
                        for (i,j),d in minors.items())
            square += times[p]*times[q]*minor**2
        voltage = sum(k*prod(segments[i] for i in ex) for ex,k in self.voltage_terms)
        poly = square*voltage
        ordered = 0
        for ex,c in poly.terms():
            denominator = 1
            for group in self.groups:
                degree = 0
                for i in group:
                    degree += int(ex[i])+1
                    denominator *= degree
            ordered += int(c)*pow(denominator,-1,PRIME)
        if check:
            # In independent gap lengths each ordered chamber is a standard simplex.
            coordinates = [sum(x[h] for h in group[:i+1]) for group in self.groups for i in range(len(group))]
            gaps = poly.compose(*coordinates)
            direct = 0
            for ex,c in gaps.terms():
                ex = tuple(map(int,ex))
                denominator = prod(factorial(len(group)+sum(ex[i] for i in group)) for group in self.groups)
                direct += int(c)*prod(map(factorial,ex))*pow(denominator,-1,PRIME)
            assert (direct-ordered)%PRIME == 0
            if times == METRIC:
                self.direct(times,square,voltage)
        # The voltage polynomial is four times the mean. Both soft label orders contribute.
        factor = self.weight/2*prod(times[e] for e in self.marks)
        return ordered*residue(factor)%PRIME


def star_period(color=COLOR):
    from flint import fmpq, fmpz_mpoly_ctx
    from fourfold_quartic_prism import project_squares
    from quartic_seven_plane_prism import cut_weights
    from quartic_seven_plane_prism_check import dualize, kirchhoff, pair_graph
    from quartic_seven_plane_single import cuts

    probe = PrimitiveChamber(3,(3,3,3,7,7,7),color)
    ctx = fmpz_mpoly_ctx.get(tuple(f't{i}' for i in range(15)))
    times = ctx.gens()
    paths = []
    for e in range(9):
        start = sum(map(len,paths))
        paths.append(tuple(range(start,start+probe.marks.count(e)+1)))
    lengths = [sum(times[i] for i in p) for p in paths]
    cycles = probe.cycles
    minors = {(i,j):numerator(probe.graph,probe.sheared,probe.flow,[(3,i),(3,j)])
              for i,j in combinations(range(4),2)}
    if not any(minors.values()):
        assert not probe.numerator(METRIC,check=True)
        return F(),F()
    assert all(not d for pair,d in minors.items() if pair!=(1,3))
    # The sole surviving omitted pair is (1,3) in the cycle basis of G(6,CORE).
    i,j = 1,3
    weight = sum((cycles[i][p]*cycles[j][q]-cycles[i][q]*cycles[j][p])**2*lengths[p]*lengths[q]
                 for p,q in combinations(range(9),2))
    voltage = sum(k*prod(times[i] for i in ex) for ex,k in probe.voltage_terms)
    data = (probe.graph,ctx,None,lengths,None,paths,None)
    transformed = project_squares(data,[(minors[i,j],weight)],voltage,1)
    if color==COLOR:
        assert len(transformed) == 140278
    for metric in (METRIC,(11,5,2,7,9,6,8,3,4)):
        original = probe.numerator(metric,check=True)
        new,jac = dualize(metric)
        new_u = kirchhoff(pair_graph(),new)
        new_n = F(str(transformed(*(fmpq(str(v)) for v in new))))
        assert residue(probe.weight/2*new_n/new_u**13) == (
            original * pow(probe.kirchhoff(metric),-13,PRIME) * residue(jac))%PRIME
    bank = cut_weights(transformed)
    if color==COLOR:
        assert len(bank)==741
    a,b = (v*probe.weight/2 for v in cuts(bank))
    if color==COLOR:
        assert (a,b)==(-F(130303752962524429,449417733427049587200000000),F(39016635509,161759973158784000000))
        partial = sum(F(1,j**3) for j in range(1,513))
        lo,hi = partial+F(1,2*513**2),partial+F(1,2*512**2)
        assert -197 < (a+b*lo)*10**17 < (a+b*hi)*10**17 < -196
    return a,b


def shifts():
    from quartic_bipartite_shift import graph_check
    from quartic_prism_shift import symbols, taylor_derivative
    table = json.loads((Path(__file__).resolve().parent.parent/'quartic_bipartite_shift.json').read_text())
    def at_strength(row):
        return sum(c*(-12)**j for j,c in enumerate(row))
    den = at_strength(table['common_denominator'])
    assert den == -5592001781747700000
    terms = {(tuple(row['matrix']),tuple(row['derivatives'])):F(at_strength(row['numerator']),den)
             for row in table['terms']}
    graph_check(-12,terms)
    u,_,_,d,h,k = symbols()
    for times in ([1]*9,METRIC):
        lhs = F(str(u(*times)))*taylor_derivative(u,h*k,-12,times)
        rhs = 16698*taylor_derivative(u,d,-13,times)
        assert lhs==rhs


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--prism-period',action='store_true',help='Integrate the triple/triple chamber with two empty stars')
    parser.add_argument('--shifts',action='store_true',help='Check the two existing dimension identities at this power')
    args = parser.parse_args()
    for index,reference in ((3,(190,824,159,917)),(4,(323,337,36,706))):
        probe = PrimitiveChamber(index)
        assert probe.weight==8
        root = probe.divisor()
        first = probe.numerator(METRIC,check=True)
        singular = probe.numerator(root,check=True)
        orbit = 0
        for action in probe.actions:
            times = tuple(root[e] for e,_ in action)
            assert not probe.kirchhoff(times)
            orbit += probe.numerator(times)
        orbit = orbit*pow(len(probe.actions),-1,PRIME)%PRIME
        assert (first,root[0],singular,orbit) == reference
        print('Core',index,'positive metric, root, residue, symmetry mean:',*reference,flush=True)
    if args.prism_period:
        print('Two-star prism period in the rational and zeta(3) basis:',star_period(),flush=True)
    if args.shifts:
        shifts()
