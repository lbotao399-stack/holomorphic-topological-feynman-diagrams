"""Isolated genus-four double-soft fourfold Ward projection and its genus-three pentagon."""
from fractions import Fraction as F
from functools import cache
from itertools import combinations
from math import comb, factorial
from pathlib import Path
import json
import sys

import numpy as np
from flint import fmpq, fmpq_mat, fmpz_mat, fmpz_mpoly_ctx
from fourfold_quartic_basis import ORDINARY_WORDS, minors, permuted_minors, square_orbits, square_values
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from elementary_scalar import sextics
from marked_theta import green_moments

AXES=tuple(tuple(int(i==j) for j in range(4)) for i in range(4))
V=(0,-1,1,0)
POSITIONS=fmpz_mpoly_ctx.get(('s0','s1','s2','s3')).gens()
COEFFICIENTS=json.loads((Path(__file__).resolve().parent.parent/'elementary_scalar.json').read_text())['three_collinear_sextics']


def balance(rows):
    return tuple(-sum(values) for values in zip(*rows))


def frame(t):
    soft=tuple(t*x for x in V)
    hard=[*AXES,AXES[3]]
    return [*hard,balance([*hard,soft,soft]),soft,soft]


def simplex(t,direction):
    soft=tuple(t*x for x in direction)
    return [*AXES,balance([*AXES,soft,soft]),soft,soft]


def trace_tests():
    indices,signs=permuted_minors()
    groups=square_orbits(indices)

    @cache
    def traces(rows):
        m=minors(rows)
        values=square_values(groups,m)
        signed=np.asarray(m,dtype=object)[indices]*signs
        for word in ORDINARY_WORDS:
            values.append(int(np.prod(signed[:,word],axis=1).sum()))
        return values

    result=[F(0)]*43
    for t in range(1,6):
        hard=frame(t)
        rows=traces(tuple(simplex(t,(1,1,0,0))))
        values=[16*x for x in rows]
        for i,sign in ((0,1),(1,-1)):
            for j,row in enumerate(hard):
                if not row[3] or i==j: continue
                merged=tuple(x+y for x,y in zip(hard[i],row))
                fused=tuple(p for k,p in enumerate(hard) if k not in (i,j))+(merged,)
                for k,x in enumerate(traces(fused)):
                    values[k]+=3*sign*row[3]*x
        # Every minor is affine in t and the two soft slots each contribute at least t^2.
        weight=F((-1)**(t-1)*comb(5,t),t**4)
        result=[x+weight*y for x,y in zip(result,values)]
    assert not any(result)
    print('Ward projection on all 41 square orbits and both signed trace words: -16/3 times the mixed coefficient.')


@cache
def green(powers,i,j):
    return green_moments(powers,i,j)[1]


def pentagon(rows,plus,minus,direction):
    """Choose the four singletons, ground the six-point genus-three vertex, and sum label orders."""
    integrals=[F(0)]*3
    for singles in combinations(range(8),4):
        if plus not in singles and minus not in singles:
            continue
        rays=[rows[i] for i in singles]
        determinant=int(fmpz_mat(rays).det())
        if not determinant: continue
        group=[row for i,row in enumerate(rows) if i not in singles]
        strengths=[sum(x*y for x,y in zip(p,direction)) for p in rays]
        # Six conserved vertex rows: the four grouped waves, the cycle current
        # u=sum s_i p_i, and its complementary current sum (1-s_i) p_i.
        # Their signed four-row minors are differences of six Gale coordinates.
        reference=int(fmpz_mat(group).det())
        gale=[]
        for j in range(4):
            gale.append(sum(x*int(fmpz_mat(group[:j]+[r]+group[j+1:]).det())
                            for x,r in zip(POSITIONS,rays)))
        gale.extend((-reference,0))
        squared={(i,j):(gale[i]-gale[j])**2 for i,j in combinations(range(6),2)}
        for idx,polynomial in enumerate(sextics(squared)):
            if not polynomial:continue
            for powers,coefficient in polynomial.terms():
                moment=F(0)
                powers=tuple(map(int,powers))
                for i,original in enumerate(singles):
                    sign=(original==plus)-(original==minus)
                    if sign:
                        moment+=sign*sum(x*green(powers,i,j) for j,x in enumerate(strengths))
                integrals[idx]+=F(int(coefficient)*determinant**2,2)*moment
    # fourfold_primary uses the minus sign and the six labeled Abel assignments.
    return tuple(-sum(a*F(b) for a,b in zip(integrals,column))/6 for column in zip(*COEFFICIENTS))


def mixed_pentagon():
    # The integrand is divisible by t^4 by two independent quadratic soft
    # orders unless both soft rows are singletons, when the determinant is zero.
    # A grouped soft row is quadratically soft in the genus-three action.
    readings=[pentagon(frame(t),0,1,AXES[3]) for t in range(1,6)]
    weights=[F((-1)**(t-1)*comb(5,t),t**4) for t in range(1,6)]
    value=tuple(sum(w*c for w,c in zip(weights,col)) for col in zip(*readings))
    # Independent evaluation at t=-1 by the degree-four interpolation of V(t)/t^4.
    negative=tuple(sum((-1)**(i-1)*comb(5,i)*F(6*i,i+1)*c/i**4 for i,c in enumerate(col,1))
                   for col in zip(*readings))
    assert negative==pentagon(frame(-1),0,1,AXES[3])
    assert value==(F(740466719,23832365356800000),-F(469153,18389170800000))
    return value


def bare_minor(rays,cycles,lengths,first,second):
    """The two differentiated jump columns in the unsigned projective cofactor."""
    rows=[[1]*len(lengths)+[0,0]]
    for mu in range(4):
        for k,cycle in enumerate(cycles):
            rows.append([str(c*p[mu]) for c,p in zip(cycle,rays)]
                        +[str(first[k]*V[mu]),str(second[k]*V[mu])])
    return F(str(fmpq_mat(rows).det()))/sum(lengths)


def integrated_soft_placements(rays, cycles, lengths, direction, edge_readings=None):
    """Sum the squared cofactor, optionally testing one soft mark with an affine voltage."""
    genus = len(cycles)
    remaining = len(lengths)
    removed = 1 + len(direction)*genus - remaining
    if removed < 1:
        raise ValueError("Expected a hard graph with at least one soft mark removed")
    if removed > genus:
        return F(0)
    hard = [[1]*remaining]
    jumps = [[0]*genus]
    for component, strength in enumerate(direction):
        for k, cycle in enumerate(cycles):
            hard.append([str(c*p[component]) for c,p in zip(cycle,rays)])
            jumps.append([str(strength*int(k==j)) for j in range(genus)])

    reduced, rank = fmpq_mat(hard).transpose().rref()
    if rank < remaining:
        return F(0)
    pivots = [next(j for j in range(len(hard)) if reduced[i,j])
              for i in range(remaining)]
    residual = [i for i in range(len(hard)) if i not in pivots]
    hard_pivot = fmpq_mat([hard[i] for i in pivots])
    jump_pivot = fmpq_mat([jumps[i] for i in pivots])
    projection = (fmpq_mat([jumps[i] for i in residual])
                  - fmpq_mat([hard[i] for i in residual])*hard_pivot.solve(jump_pivot))
    projected_cycles = projection*fmpq_mat(cycles)

    def weighted_gram(weights):
        rows = projected_cycles.tolist()
        return fmpq_mat([[sum(a*b*w for a,b,w in zip(c,d,weights))
                          for c in rows] for d in rows])

    gram = weighted_gram([fmpq(str(t)) for t in lengths])
    norm = hard_pivot.det()**2 * gram.det()
    if norm and edge_readings is not None:
        if len(edge_readings)!=remaining:
            raise ValueError("Expected one mean voltage for each hard segment")
        response = gram.solve(weighted_gram([fmpq(str(t*x)) for t,x in zip(lengths,edge_readings)]))
        norm *= sum(response[k,k] for k in range(removed))/removed
    return factorial(removed)*F(str(norm))/sum(lengths)**2


def hard_bipartite(marks):
    core=[(i,j) for i in range(3) for j in range(3,6)]
    edges=[]
    parents=[]
    for e,(left,right) in enumerate(core):
        chain=[left,*[6+i for i,p in enumerate(marks) if p==e],right]
        edges.extend(zip(chain,chain[1:]))
        parents.extend([e]*(len(chain)-1))
    cycles=[]
    for i in (1,2):
        for j in (4,5):
            rectangle={(0,3):1,(i,3):-1,(i,j):1,(0,j):-1}
            cycles.append([rectangle.get(core[e],0) for e in parents])
    lengths=[F(2+(i*i+3*i)%17) for i in range(len(edges))]
    return edges,cycles,lengths


def bare_check():
    from point_graph import solve_graph
    from scalar_fan import cofactor

    edges,cycles,lengths=hard_bipartite((0,1,2,3,4,8))
    zero=(0,)*4
    def sources(t):
        return [*([zero]*6),*frame(t)]
    rays,potentials=solve_graph(edges,sources(0)[:12],lengths,0)
    mean_voltage = [(potentials[i][3]+potentials[j][3])/2 for i,j in edges]
    cols=list(zip(*cycles))
    skew=fmpq_mat(4,4)
    for i,j in combinations(range(4),2):
        value=bare_minor(rays,cycles,lengths,AXES[i],AXES[j])
        skew[i,j],skew[j,i]=str(value),str(-value)
    gram=fmpq_mat([[str(sum(a*b*t for a,b,t in zip(r,s,lengths))) for s in cycles] for r in cycles])
    km=skew*gram
    all_pairs=F(0)
    soft_voltage=F(0)
    picked=None
    for i,j in combinations(range(15),2):
        lead=bare_minor(rays,cycles,lengths,cols[i],cols[j])
        assert lead==sum(F(str(skew[k,l]))*cols[i][k]*cols[j][l] for k in range(4) for l in range(4))
        pair_weight=2*lengths[i]*lengths[j]*lead**2
        all_pairs+=pair_weight
        soft_voltage+=pair_weight*(mean_voltage[i]+mean_voltage[j])/2
        if lead and picked is None:picked=(i,j,lead)
    square=km*km
    assert all_pairs==-sum(F(str(square[i,i])) for i in range(4)) and all_pairs>0
    assert integrated_soft_placements(rays,cycles,lengths,V)==all_pairs
    assert soft_voltage and integrated_soft_placements(rays,cycles,lengths,V,mean_voltage)==soft_voltage
    assert skew[0,1]*skew[2,3]-skew[0,2]*skew[1,3]+skew[0,3]*skew[1,2]==0
    assert skew.rank()==2
    first,second,lead=picked
    assert lead
    # Check two different positions, restoring both labels and the balanced recoil.
    for a,b in ((F(1,3),F(2,5)),(F(4,7),F(7,8))):
        full_edges,times,parents=[],[],[]
        for i,(l,r) in enumerate(edges):
            if i==first or i==second:
                v,x=(12,a) if i==first else (13,b)
                full_edges.extend(((l,v),(v,r)));times.extend((lengths[i]*x,lengths[i]*(1-x)))
                parents.extend((i,i))
            else:
                full_edges.append((l,r));times.append(lengths[i]);parents.append(i)
        lifted=[[row[i] for i in parents] for row in cycles]
        # Keep the signed cofactor here. It has degree <=4 since every t-dependent
        # current is in the single target direction v and there are four cycles.
        def coefficient(t):
            y,restored_potential=solve_graph(full_edges,sources(t),times,0)
            for node,edge,fraction in ((12,first,a),(13,second,b)):
                left,right=edges[edge]
                assert restored_potential[node][3]==(1-fraction)*potentials[left][3]+fraction*potentials[right][3]
            rows=[[1]*17]+[[str(c*x[mu]) for c,x in zip(cycle,y)] for mu in range(4) for cycle in lifted]
            signed=F(str(fmpq_mat(rows).det()))/sum(times)
            assert abs(signed)==cofactor(y,lifted,times)
            return signed
        linear=coefficient(1)-coefficient(-1)
        quadratic=(16*(coefficient(1)+coefficient(-1))-coefficient(2)-coefficient(-2))/24
        assert quadratic**2==lead**2
        assert 8*linear==coefficient(2)-coefficient(-2)
    print('Two-soft bare minor on two sets of rational insertion positions:',first,second,lead**2)


def check_more_soft_slots():
    from math import prod
    from point_graph import solve_graph
    from scalar_fan import cofactor

    for marks in ((0,1,3,4,8),(0,1,3,4)):
        edges,cycles,lengths = hard_bipartite(marks)
        hard = AXES[:len(marks)-1]
        flows,potential = solve_graph(edges,[(0,)*4]*6+[*hard,balance(hard)],lengths,0)
        count = 8-len(marks)
        value = integrated_soft_placements(flows,cycles,lengths,AXES[3])
        assert value > 0
        if count==4:
            gram = fmpq_mat([[str(sum(a*b*t for a,b,t in zip(c,d,lengths)))
                              for c in cycles] for d in cycles])
            reduced = cofactor([row[:3] for row in flows],cycles,lengths)
            assert value == factorial(count)*F(str(gram.det()))*reduced**2
            strengths = (2,-1,3,5)
            shifts = (1,-2,4,-3-sum(strengths))
            source = [(0,)*4]*6+[(*p[:3],q) for p,q in zip([*hard,balance(hard)],shifts)]
            source.extend((0,0,0,a) for a in strengths)
            assert not any(map(sum,zip(*source)))
            cotrees = [(s,fmpz_mat([[c[i] for i in s] for c in cycles]).det())
                       for s in combinations(range(len(edges)),count)]
            supported = [(s,d) for s,d in cotrees if d]
            for chosen,jump in (supported[0],supported[-1]):
                full_edges,times,parents = [],[],[]
                for i,(left,right) in enumerate(edges):
                    if i in chosen:
                        j = chosen.index(i)
                        position = F(j+1,j+3)
                        full_edges.extend(((left,10+j),(10+j,right)))
                        times.extend((lengths[i]*position,lengths[i]*(1-position)))
                        parents.extend((i,i))
                    else:
                        full_edges.append((left,right))
                        times.append(lengths[i])
                        parents.append(i)
                full_flow,full_potential = solve_graph(full_edges,source,times,0)
                actual = cofactor(full_flow,[[c[i] for i in parents] for c in cycles],times)
                assert actual == abs(int(jump)*prod(strengths))*reduced
                assert all(full_potential[i][:3]==v[:3] for i,v in potential.items())
        else:
            # Literal augmented cofactor for every set of three target jumps.
            terms = F(0)
            for selected in combinations(range(len(edges)),count):
                rows = [[1]*len(edges)+[0]*count]
                for mu in range(4):
                    for c in cycles:
                        rows.append([str(k*p[mu]) for k,p in zip(c,flows)]
                                    +[int(mu==3)*c[i] for i in selected])
                minor = F(str(fmpq_mat(rows).det()))/sum(lengths)
                terms += prod(lengths[i] for i in selected)*minor**2
            assert value == factorial(count)*terms
    print('Three soft slots: all 364 triples. Four soft slots: exact threefold reduction.')


def leading_pentagon():
    """One soft at the bulk vertex, one at a singleton; extract t^4 before integration."""
    hard=frame(0)[:6]
    positions=fmpz_mpoly_ctx.get(('s0','s1','s2')).gens()
    total=F(0)
    for singles in combinations(range(6),3):
        if 0 not in singles and 1 not in singles:continue
        rays=[hard[i] for i in singles]
        volume=int(fmpz_mat([*rays,V]).det())
        if not volume: continue
        group=[q for i,q in enumerate(hard) if i not in singles]

        def with_current(fixed):
            return sum(s*int(fmpz_mat([*fixed,q]).det()) for s,q in zip(positions,rays))

        d=with_current(group)
        replacement=[with_current(group[:i]+[V]+group[i+1:]) for i in range(3)]
        replacement.append(int(fmpz_mat([*group,V]).det()))
        # Cramer's rule for the soft wave on the hard simplex. The fifth affine
        # coefficient is zero. The two choices of the singleton soft label
        # cancel the cycle-reflection factor 1/2.
        symbol=d**4*(5*sum(x*x for x in replacement)-sum(replacement)**2)
        for powers,coefficient in symbol.terms():
            powers=tuple(map(int,powers))
            for i,label in enumerate(singles):
                sign=(label==0)-(label==1)
                if sign:
                    moment=sum(q[3]*green(powers,i,j) for j,q in enumerate(rays))
                    total+=F(sign*volume**2*int(coefficient),4)*moment
    assert total==F(4,225)
    return total


def check_pentagon():
    normal=((1,0,0),(0,1,0),(0,0,1),(-1,0,1),(0,-1,-2))
    old=[(a,0,0,0) for a in (1,1,-2)]+[(0,*v) for v in normal]
    old_value=pentagon(old,3,4,AXES[3])
    assert old_value==(F(3,280)*F(22599229,252194342400),-F(3,280)*F(143249,1945944000))
    print('Reference zero-sum triple-collinear pentagon reproduced.')
    measured=mixed_pentagon()
    axial_three=(F(740466719,7565830272000*56),-F(469153,5837832000*56))
    assert measured==tuple(leading_pentagon()*c for c in axial_three)
    print('Leading three-position pentagon:',leading_pentagon(),'times the genus-three axial soft coefficient.')
    print('Mixed two-soft pentagon in 1,zeta(3) units:',measured)


if __name__=='__main__':
    trace_tests()
    check_pentagon()
    bare_check()
    check_more_soft_slots()
