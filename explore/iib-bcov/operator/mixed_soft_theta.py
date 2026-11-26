"""The genus-two vertex on the theta graph in the mixed genus-four soft Ward row."""
from fractions import Fraction
from itertools import combinations
from math import prod

from quartic_seven_plane_theta import (CTX, DEN, S, W, chambers, cofactor,
                                      integrate_chamber)
from flint import fmpq_mat, fmpz_mat


HARD = ((1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1),(0,0,0,1),(-1,-1,-1,-2))


def projected(p):
    # The removed soft direction is e3-e2; these three coordinates annihilate it.
    return p[0], p[1]+p[2], p[3]


def contract_current(mark, p, path, chamber, kind):
    """Four-volume at the marked vertex with all its hard inputs held fixed."""
    if kind == 'interior':
        minors = [int(fmpz_mat([*mark,q]).det()) for q in p]
        group = chamber[path[0]]
        place = group.index(0)
        value = W[path[0]]*sum(s*x for s,x in zip(S,minors))
        return value + sum((int(place<j)-S[i])*minors[i] for j,i in enumerate(group))

    minors = [[int(fmpz_mat([*mark,pi,pj]).det()) for pj in p] for pi in p]
    h0,h1 = chamber[:2]
    base = sum(((1-S[i])*(1-S[j])*minors[i][j] for i in h0 for j in h1), CTX.constant(0))
    left = sum((S[i]*(1-S[j])*minors[i][j] for i in range(4) for j in h1), CTX.constant(0))
    right = sum(((1-S[i])*S[j]*minors[i][j] for i in h0 for j in range(4)), CTX.constant(0))
    return base + W[0]*left + W[1]*right


def group_integral(hard, target, marked, kind):
    others = [i for i in range(6) if i not in marked]
    roots = [(i,) for i in others]
    if kind=='interior':
        roots.insert(0,marked)
    p = [tuple(map(sum,zip(*(hard[i] for i in root)))) for root in roots]
    charges = [q[3] for q in p]
    signs = [int(target[0] in root)-int(target[1] in root) for root in roots]
    mark = [hard[i] for i in marked]
    if not any(signs) or fmpz_mat(mark).rank()<len(mark):
        return Fraction(0)
    normal = [projected(q) for q in p]
    if fmpz_mat(normal).rank()<3:
        return Fraction(0)
    global_voltage = sum(s*x for s,x in zip(S,charges))*sum(s*x for s,x in zip(S,signs))
    result = 0
    for path,order in chambers():
        det = cofactor(normal,path)
        if not det:
            continue
        vertex = contract_current(mark,p,path,order,kind)
        if not vertex:
            continue
        density = det**2*vertex**4
        result += integrate_chamber(density*global_voltage,order)
        for arm,group in enumerate(order):
            local_voltage = sum((signs[i]*charges[j]*(S[group[min(ii,jj)]]-S[i]*S[j])
                                 for ii,i in enumerate(group) for jj,j in enumerate(group)), CTX.constant(0))
            if local_voltage:
                result += integrate_chamber(density*local_voltage,order,arm)
    return Fraction(result,DEN)


def theta_component(hard, target, kind):
    count = 2 if kind=='branch' else 3
    total = sum(group_integral(hard,target,marks,kind) for marks in combinations(range(6),count))
    # Two labeled soft slots give 2!, the five-field vertex is det^4/62208,
    # and an interior distinguished vertex admits interchange of the theta endpoints.
    return total * Fraction(2, 62208*(1 if kind=='branch' else 2))


def check_restored_currents(kind):
    """Compare the soft placement sum to literal nine-edge, four-dimensional networks."""
    from point_graph import solve_graph
    from quartic_seven_plane_theta import evaluate
    from scalar_fan import cofactor as current_cofactor

    marked = (0,2) if kind=='branch' else (0,2,3)
    roots = [(i,) for i in range(6) if i not in marked]
    if kind=='interior':
        roots.insert(0,marked)
    mark = [HARD[i] for i in marked]
    event = [tuple(map(sum,zip(*(HARD[i] for i in group)))) for group in roots]
    zero = (0,)*4
    source = [tuple(map(sum,zip(*mark))) if kind=='branch' else zero,zero,*event]
    vertex = 0 if kind=='branch' else 2
    sites = {i:j+2 for j,group in enumerate(roots) for i in group}
    if kind=='branch':
        sites.update(dict.fromkeys(marked,0))
    path = (0,0,1,1)
    chamber = ((0,1),(2,3),())
    s = [Fraction(k,13) for k in (2,5,3,6)]
    w = [Fraction(k,10) for k in (2,3,5)]
    edges, lengths, arms = [], [], []
    for h,group in enumerate(chamber):
        sites_on_arm = [0,*[2+i for i in group],1]
        positions = [0,*[s[i] for i in group],1]
        edges.extend(zip(sites_on_arm,sites_on_arm[1:]))
        lengths.extend((b-a)/w[h] for a,b in zip(positions,positions[1:]))
        arms.extend([h]*(len(group)+1))
    cycles = [[int(h==k)-int(h==1) for h in arms] for k in (0,2)]
    tree = sum(1/(w[i]*w[j]) for i,j in combinations(range(3),2))
    currents, potentials = solve_graph(edges,source,lengths,0)
    inside = [[(int(a==vertex)-int(b==vertex))*x for x in row]
              for (a,b),row in zip(edges,currents) if vertex in (a,b)]
    volume = fmpq_mat([[str(v) for v in q] for q in [*mark,*inside[:4-len(mark)]]]).det()
    assert Fraction(str(volume))**4==evaluate(contract_current(mark,event,path,chamber,kind),s+w)**4
    voltage = potentials[sites[0]][3]-potentials[sites[1]][3]
    expected = Fraction()
    for i,group in enumerate(roots):
        sign = int(0 in group)-int(1 in group)
        for j,q in enumerate(event):
            green = s[i]*s[j]
            if path[i]==path[j]:
                green += (min(s[i],s[j])-s[i]*s[j])/w[path[i]]
            expected += sign*q[3]*green
    assert voltage==expected and volume and voltage
    hard_cofactor = current_cofactor([projected(q) for q in currents],cycles,lengths)
    assert hard_cofactor==abs(evaluate(cofactor([projected(q) for q in event],path),s+w))/(prod(w)*tree**2)
    soft = (0,-1,1,0)
    placement_sum = 0
    for i,j in combinations(range(len(edges)),2):
        jump = cycles[0][i]*cycles[1][j]-cycles[0][j]*cycles[1][i]
        for t,left,right in ((1,Fraction(1,3),Fraction(2,5)),(-2,Fraction(3,5),Fraction(1,4))):
            restored, times, parent = [], [], []
            for k,(a,b) in enumerate(edges):
                if k in (i,j):
                    cut, site = (left,6) if k==i else (right,7)
                    restored.extend(((a,site),(site,b)))
                    times.extend((cut*lengths[k],(1-cut)*lengths[k]))
                    parent.extend((k,k))
                else:
                    restored.append((a,b))
                    times.append(lengths[k])
                    parent.append(k)
            waves = source+[tuple(t*x for x in soft)]*2
            waves[sites[5]] = tuple(x-2*t*v for x,v in zip(waves[sites[5]],soft))
            flow, phi = solve_graph(restored,waves,times,0)
            observed = current_cofactor(flow,[[c[k] for k in parent] for c in cycles],times)
            assert observed==t**2*abs(jump)*hard_cofactor
            assert phi[sites[0]][3]-phi[sites[1]][3]==voltage
        placement_sum += 2*lengths[i]*lengths[j]*(jump*hard_cofactor)**2
    assert placement_sum==2*tree*hard_cofactor**2 and placement_sum
    return voltage, Fraction(str(volume))**4


def check_three_line_normalization():
    """Recover the branch and interior readings on the three-collinear fourfold row."""
    # The fixed three-line fixture in fourfold_quartic_ward.lower_vertices().
    normal = ((1,0,0),(0,1,0),(0,0,1),(-1,0,1),(0,-1,-2))
    strengths = (1,1,-2)
    soft = (0,-1,1,0)
    hard = [soft,*[(x,y,0,z) for x,y,z in normal]]
    hard[1] = tuple(x-y for x,y in zip(hard[1],soft))
    # The isolated unit hard line must enter the action vertex. Its projection is zero.
    # With all three line slots restored, the chosen ordinary pair enters quadratically
    # and the remaining line enters the fourth power of the action determinant.
    line_factor = prod(strengths)**2*sum(c*c for c in strengths)
    assert line_factor==24
    reduced = (Fraction(-159293,91683900000),Fraction(-1423,63504000))
    for kind,coefficient in zip(('branch','interior'),reduced):
        assert line_factor*theta_component(hard,(1,2),kind)==coefficient/1296


if __name__ == '__main__':
    for kind,answer in (('branch',Fraction(896261,3802314700800000)),
                        ('interior',Fraction(1,23702740992))):
        print(kind,"network",check_restored_currents(kind),flush=True)
        observed = theta_component(HARD,(0,1),kind)
        assert observed==answer
        print(kind,observed,flush=True)
    check_three_line_normalization()
