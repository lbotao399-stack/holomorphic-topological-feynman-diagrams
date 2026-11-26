"""Integrate all six freely placed hard waves and both soft waves at fixed core lengths."""
from fractions import Fraction as F
from functools import cache
from itertools import combinations, permutations, product
from math import factorial, prod

from flint import fmpq, fmpq_mat, fmpq_mpoly_ctx, fmpq_poly
from fourfold_mixed_soft_ward import V, hard_bipartite, integrated_soft_placements
from mixed_bare_ranks import WAVES
from point_graph import solve_graph
from quartic_seven_plane_double import polynomial_det


CTX=fmpq_mpoly_ctx.get(tuple(f'x{i}' for i in range(6)))
X=CTX.gens()
SOFT=tuple(combinations(range(4),2))


def cofactor(means,jumps):
    # Move the eight constant jump columns in front of the nine mean columns.
    transposed=fmpq_mat(jumps).transpose()
    echelon,rank=transposed.rref()
    if rank<8:return CTX.constant(0)
    selected=[next(j for j in range(17) if echelon[i,j]) for i in range(8)]
    others=[i for i in range(17) if i not in selected]
    pivot=fmpq_mat([jumps[k] for k in selected])
    cancel=fmpq_mat([jumps[k] for k in others])*pivot.inv()
    remaining=[[means[i][j]-sum(cancel[k,l]*means[s][j] for l,s in enumerate(selected) if cancel[k,l])
                for j in range(9)] for k,i in enumerate(others)]
    assert others[0]==0 and all(v==1 for v in remaining[0])
    # Subtract the first mean column from each of the other eight.
    reduced=[[v-row[0] for v in row[1:]] for row in remaining[1:]]
    return (-1)**sum(s-i for i,s in enumerate(selected))*pivot.det()*polynomial_det(reduced,CTX)


def integrate_unit(poly):
    return sum(c/fmpq(prod(int(j)+1 for j in ex)) for ex,c in poly.terms())


@cache
def grounded_response(edges,times):
    nodes=1+max(max(edge) for edge in edges)
    identity=[[int(i==j) for j in range(nodes)] for i in range(nodes)]
    return solve_graph(edges,identity,times,0)


def voltage_weights(averaged):
    if not averaged:
        return [(i,j,F(sign*wave[3])) for i,sign in ((0,1),(1,-1)) for j,wave in enumerate(WAVES[:6]) if wave[3]]
    weights=[]
    for i,j in permutations((0,5)):
        weights.append((i,j,F(-1)))
    for i,j in product((0,5),(3,4)):weights.append((i,j,F(1,2)))
    for i,j in product((1,2),(0,5)):weights.append((i,j,F(1,2)))
    for i,j in product((1,2),(3,4)):weights.append((i,j,F(-1,2)))
    return weights


class HardCube:
    def __init__(self,edges,cycles,times,marks,averaged=False):
        self.edges,self.cycles,self.times,self.marks=edges,cycles,times,marks
        self.waves=WAVES[:6]
        self.jacobian=prod(times[e] for e in marks)
        currents,potentials=grounded_response(tuple(edges),tuple(times))
        def mean(row,mu):
            return sum(fmpq(str(row[u]))*q[mu]+fmpq(str(row[v]-row[u]))*q[mu]*x
                       for (u,v),q,x in zip((edges[e] for e in marks),self.waves,X) if q[mu])
        flow=[[mean(row,mu) for mu in range(4)] for row in currents]
        means=[[CTX.constant(1)]*9]+[[r[e]*flow[e][mu] for e in range(9)] for mu in range(4) for r in cycles]
        js=[[0]*6]+[[r[e]*q[mu] for e,q in zip(marks,self.waves)]
                       for mu in range(4) for r in cycles]
        minors=[cofactor(means,[[*row,0,0] if k==0 else
                   [*row,V[(k-1)//4]*int((k-1)%4==a),V[(k-1)//4]*int((k-1)%4==b)]
                   for k,row in enumerate(js)]) for a,b in SOFT]
        gram=[[sum(t*i*j for i,j,t in zip(u,v,times)) for v in cycles] for u in cycles]
        self.square=2*sum((p*q*fmpq(str(gram[a][c]*gram[b][d]-gram[a][d]*gram[b][c]))
                    for (a,b),p in zip(SOFT,minors) if p for (c,d),q in zip(SOFT,minors) if q),
                    CTX.constant(0))/fmpq(str(sum(times)))**2
        self.voltage=CTX.constant(0)
        self.green=[]
        for i,j,weight in voltage_weights(averaged):
            u,v=edges[marks[i]]
            y,z=edges[marks[j]]
            read_u=(1-X[j])*fmpq(str(potentials[u][y]))+X[j]*fmpq(str(potentials[u][z]))
            read_v=(1-X[j])*fmpq(str(potentials[v][y]))+X[j]*fmpq(str(potentials[v][z]))
            self.voltage+=fmpq(str(weight))*((1-X[i])*read_u+X[i]*read_v)
            if marks[i]==marks[j]:
                self.green.append((i,j,weight*times[marks[i]]))

    def cube(self):
        result=integrate_unit(self.square*self.voltage)
        for ex,c in self.square.terms():
            powers=list(map(int,ex))
            weight=c/fmpq(prod(i+1 for i in powers))
            for i,j,scale in self.green:
                a,b=powers[i],powers[j]
                if i==j:moment=fmpq(a+1,(a+2)*(a+3))
                else:moment=fmpq((a+1)*(b+1),(a+2)*(b+2)*(a+b+3))
                result+=weight*fmpq(str(scale))*moment
        return F(str(result))*self.jacobian

    def ordinary(self,points):
        edges=[];lengths=[];parents=[]
        for e,((u,v),time) in enumerate(zip(self.edges,self.times)):
            order=sorted((points[i],6+i) for i,h in enumerate(self.marks) if h==e)
            path=[(F(),u),*order,(F(1),v)]
            for (s,i),(t,j) in zip(path,path[1:]):
                assert t>s
                edges.append((i,j));lengths.append((t-s)*time);parents.append(e)
        source=[(0,)*4]*6+list(self.waves)
        flow,phi=solve_graph(edges,source,lengths,0)
        cycles=[[r[e] for e in parents] for r in self.cycles]
        square=integrated_soft_placements(flow,cycles,lengths,V)
        potential=phi[6][3]-phi[7][3]
        qs=tuple(fmpq(str(p)) for p in points)
        assert self.square(*qs)==fmpq(str(square))
        voltage=self.voltage(*qs)+sum(fmpq(str(t))*(min(qs[i],qs[j])-qs[i]*qs[j]) for i,j,t in self.green)
        assert voltage==fmpq(str(potential))

    def ordered(self):
        # An independent integration measure: recover every ordered simplex and
        # choose the earlier of the two marks on each path before integrating.
        groups=[[i for i,k in enumerate(self.marks) if e==k] for e in sorted(set(self.marks))]
        out=fmpq(0)
        for permuted in product(*(permutations(group) for group in groups)):
            rank={i:k for group in permuted for k,i in enumerate(group)}
            voltage=self.voltage
            for i,j,scale in self.green:
                earlier=i if rank[i]<=rank[j] else j
                voltage+=fmpq(str(scale))*(X[earlier]-X[i]*X[j])
            for ex,c in (self.square*voltage).terms():
                denom=1
                for group in permuted:
                    degree=0
                    for i in group:
                        degree+=int(ex[i])+1
                        denom*=degree
                out+=c/fmpq(denom)
        return F(str(out))*self.jacobian


def labeled_cubes(probe,marks,times,averaged=False):
    from collections import Counter

    stabilizer=[g for g in probe.actions if tuple(sorted(g[e][0] for e in marks))==marks]
    isometries=[g for g in stabilizer if all(times[g[e][0]]==t for e,t in enumerate(times))]
    counts=Counter()
    for labeled in set(permutations(marks)):
        images=[tuple(g[e][0] for e in labeled) for g in isometries]
        if averaged:
            def ordered_pairs(r):
                u=sorted((r[0],r[5]))
                t=sorted((r[1],r[2]))
                z=sorted((r[3],r[4]))
                return u[0],*t,*z,u[1]
            counts[min(ordered_pairs(r) for r in images)]+=1
        else:
            # The two fourth-coordinate waves have identical source and observation.
            counts[min(images+[(*r[:3],r[4],r[3],r[5]) for r in images])]+=1
    return counts,len(stabilizer)


def placement(probe,counts,stabilizer,times,averaged=False):
    full=F()
    for representative,count in counts.items():
        cube=HardCube(probe.core,probe.cycles,times,representative,averaged)
        full+=count*cube.cube()
    return full/stabilizer


def ordered_numerator_check(probe,marks,times,density,core_index=4):
    from mixed_bare_primitive_pole import PRIME, PrimitiveChamber, residue
    from mixed_prism_one_star import chamber_rows

    stabilizer=[g for g in probe.actions if tuple(sorted(g[e][0] for e in marks))==marks]
    equivalent={tuple(times[g[e][0]] for e in range(9)) for g in stabilizer}
    expected=0
    for color,_ in chamber_rows(marks,core_index):
        chamber=PrimitiveChamber(core_index,marks,color)
        expected+=sum(chamber.numerator(tuple(metric)) for metric in equivalent)
    normalized=residue(density)*pow(probe.kirchhoff(times),9,PRIME)
    assert (normalized*len(equivalent)-expected)%PRIME==0


def matching_check(marks):
    from mixed_bare_primitive_pole import PrimitiveChamber

    probe=PrimitiveChamber(4,marks)
    unit=(1,)*9
    counts,stabilizer=labeled_cubes(probe,marks,unit)
    full=placement(probe,counts,stabilizer,unit)
    averaged,order=labeled_cubes(probe,marks,unit,True)
    assert full==placement(probe,averaged,order,unit,True)
    ordered_numerator_check(probe,marks,unit,full)
    print('Unit metric, placement, labeled assignments, cube orbits, signed integral:',marks,
          sum(counts.values()),len(counts),full,flush=True)


def all_unit(core_index,destination):
    from collections import Counter
    from pathlib import Path
    from mixed_bare_primitive_pole import PrimitiveChamber
    from mixed_bare_symbolic import classified_orbits
    import json

    path=Path(destination)
    records=json.loads(path.read_text()) if path.exists() else {}
    core=PrimitiveChamber(core_index)
    unit=(1,)*9
    census=Counter(r[0] for r in classified_orbits(core_index)[0])
    for marks,count in sorted(census.items()):
        word=''.join(map(str,marks))
        if word in records:continue
        assignments,stab=labeled_cubes(core,marks,unit,True)
        value=placement(core,assignments,stab,unit,True)
        records[word]=dict(ordered_color_orbits=count,unordered_assignments=len(assignments),
                           density=str(value))
        temp=path.with_suffix('.partial')
        temp.write_text(json.dumps(records,sort_keys=True,indent=2)+'\n')
        temp.replace(path)
        print(word,len(assignments),value,flush=True)
    print('Complete fixed-unit core',core_index,'placements',len(records),
          'ordered colors',sum(r['ordered_color_orbits'] for r in records.values()),
          'density',sum(F(r['density']) for r in records.values()),flush=True)


def check_all_unit(core_index,source):
    from collections import Counter
    from pathlib import Path
    from mixed_bare_primitive_pole import PrimitiveChamber
    from mixed_bare_symbolic import classified_orbits
    import json

    records=json.loads(Path(source).read_text())
    census=Counter(row[0] for row in classified_orbits(core_index)[0])
    assert set(records)=={''.join(map(str,m)) for m in census}
    unit=(1,)*9
    probe=PrimitiveChamber(core_index)
    for marks,count in sorted(census.items()):
        word=''.join(map(str,marks))
        assert records[word]['ordered_color_orbits']==count
        ordered_numerator_check(probe,marks,unit,F(records[word]['density']),core_index)
        print('Grounded chamber comparison:',word,flush=True)


def all_core_density(folder,destination):
    from collections import Counter
    from pathlib import Path
    from mixed_bare_primitive_pole import PrimitiveChamber
    from mixed_bare_symbolic import classified_orbits
    import json

    folder=Path(folder)
    unit=(1,)*9
    totals=[]
    for index in range(5):
        census=Counter(row[0] for row in classified_orbits(index)[0])
        if not census:rows={}
        elif index==4:
            rows=json.loads(Path(__file__).with_name('mixed_bipartite_unit_metric.json').read_text())
        else:
            rows=json.loads((folder/f'unit-core{index}.json').read_text())
        assert {''.join(map(str,m)):n for m,n in census.items()}=={
            word:row['ordered_color_orbits'] for word,row in rows.items()}
        u=PrimitiveChamber(index).kirchhoff(unit)
        raw=sum((F(row['density']) for row in rows.values()),F())
        totals.append(dict(core=index,kirchhoff=u,placements=len(rows),colors=sum(census.values()),
                           fixed_length_integral=str(raw),projective_density=str(raw/u**4)))
    result=dict(unit_metric=list(unit),cores=totals,placements=sum(row['placements'] for row in totals),
                colors=sum(row['colors'] for row in totals),
                projective_density=str(sum(F(row['projective_density']) for row in totals)))
    assert (result['placements'],result['colors'])==(377,19099)
    output=Path(destination)
    pending=output.with_suffix('.partial')
    pending.write_text(json.dumps(result,indent=2)+'\n')
    pending.replace(output)
    print('All cubic cores at unit lengths',result['projective_density'],flush=True)


def edge_period(marks,edge,background=(1,)*9):
    from mixed_bare_primitive_pole import PrimitiveChamber

    probe=PrimitiveChamber(4,marks)
    def times(n):return tuple(n if e==edge else t for e,t in enumerate(background))
    fixed_length=sum(times(0))
    counts,stabilizer=labeled_cubes(probe,marks,times(max(background)+1))
    def tree(weights):
        return int(fmpq_mat([[sum(x*y*t for x,y,t in zip(r,s,weights))
                              for s in probe.cycles] for r in probe.cycles]).det())
    a=tree(times(0))
    b=tree(times(1))-a
    # Current means have one affine Kirchhoff denominator. Their eight-by-eight
    # determinant enters twice. The voltage supplies one more denominator.
    # The compound soft Gram is affine in the single varying metric edge.
    bound=19+marks.count(edge)
    samples=[]
    for n in range(1,bound+3):
        density=placement(probe,counts,stabilizer,times(n))
        if n==2:ordered_numerator_check(probe,marks,times(n),density)
        value=density*(a+b*n)**17*(fixed_length+n)**2
        samples.append(fmpq(str(value)))
        print('Edge sample',n,'of',bound+2,flush=True)
    polynomial=fmpq_poly()
    binomial=fmpq_poly([1])
    for k in range(bound+1):
        polynomial+=samples[0]*binomial
        samples=[y-x for x,y in zip(samples,samples[1:])]
        binomial*=fmpq_poly([-k-1,1])/fmpq(k+1)
    assert samples==[0]
    numerator,remainder=divmod(polynomial,fmpq_poly([a,b])**8*fmpq_poly([fixed_length,1])**2)
    assert not remainder and numerator.degree()<=11
    value=sum(F(str(c))*F(a)**(k-12)*F(b)**(-k-1)*F(factorial(k)*factorial(11-k),factorial(12))
              for k,c in enumerate(numerator.coeffs()))
    if marks==(0,0,5,5,8,8) and edge==1 and background==(1,)*9:
        from sympy import Poly, Rational, Symbol, oo
        x=Symbol('x')
        p=Poly.from_list([Rational(str(c)) for c in reversed(numerator.coeffs())],x)
        assert p.count_roots(0,oo)==2
        assert p.count_roots(0,Rational(1,4))==p.count_roots(Rational(3,2),2)==1
        assert value==F(8626038015317273,361577770539338787840000000000000000)
    print('Placement, variable edge, Kirchhoff, numerator, one-edge integral:',marks,edge,(a,b),
          numerator,value,flush=True)
    return numerator,value


def check():
    edges,cycles,times=hard_bipartite(())
    for marks in ((0,0,4,4,8,8),(0,4,8,0,4,8),(0,1,2,0,4,1),(0,0,0,0,0,0)):
        cube=HardCube(edges,cycles,times,marks)
        for denominators in ((11,7,9,8,12,13),(13,12,8,9,7,11)):
            cube.ordinary(tuple(F(i+1,j) for i,j in enumerate(denominators)))
        val=cube.cube()
        if cube.square:assert cube.ordered()==val
        print('Labeled paths, cofactor-square monomials, signed cube integral:',marks,len(cube.square),val,flush=True)


if __name__=='__main__':
    import sys
    check()
    for word in sys.argv[1:]:
        if word.startswith('summary='):
            folder,output=word[8:].split(',',1)
            all_core_density(folder,output)
        elif word.startswith('check='):
            core_index,path=word[6:].split(',',1)
            check_all_unit(int(core_index),path)
        elif word.startswith('all='):
            core_index,path=word[4:].split(',',1)
            all_unit(int(core_index),path)
        elif '/' in word:
            sites,e,*background=word.split('/')
            edge_period(tuple(map(int,sites)),int(e),tuple(map(int,background[0])) if background else (1,)*9)
        else:
            matching_check(tuple(map(int,word)))
