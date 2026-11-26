"""Sum signed hard-color moments before the angular reduction of a mixed prism."""
from argparse import ArgumentParser
from collections import Counter, defaultdict
from copy import copy
from fractions import Fraction as F
from functools import cache
from itertools import combinations, combinations_with_replacement, groupby, permutations, product
from math import factorial, gcd, lcm, prod
from multiprocessing import get_context
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import os

from flint import fmpq
from fourfold_quartic_prism import CTX,T,split_paths
from mixed_prism_empty_stars import prism_data,square_roots
from mixed_prism_one_star import (
    METRIC,PRIME,PrimitiveChamber,capture,chamber_rows,image,kirchhoff,reduction,residue,
    transformed,wheel_bank,wheel_values,
)
from mixed_prism_opposite_face import bounds
from ordered_path_moments import PositionMoments
from quartic_seven_plane_double import tree_flow

HERE=Path(__file__).resolve().parent
ORIGINAL_STAR=(1,3,5)


def sparse_flow(graph,color,sheared):
    """Exchange the hard U and T pairs before choosing a particular tree current."""
    best=None
    for flip_u,flip_t in ((False,False),(True,False),(False,True),(True,True)):
        source=list(sheared)
        for group,exchange in ((0,flip_u),(1,flip_t)):
            if exchange:
                i,j=[6+k for k,c in enumerate(color) if c==group]
                source[i],source[j]=source[j],source[i]
        flow=tree_flow(graph,source)
        key=sum(bool(c) for row in flow for c in row),sum(abs(c) for row in flow for c in row)
        if best is None or key<best[0]:
            best=key,source,flow
    return best


class ColorMoments:
    def __init__(self,marks,star,core_index=3,adapted=None):
        self.adapted=os.environ.get('MIXED_ADAPTED_CYCLES')=='1' if adapted is None else adapted
        self.core_index=core_index
        self.data,self.soft=prism_data(marks,core_index)
        self.marks=marks
        self.star=star
        self.marked=tuple(i for i,path in enumerate(self.data[5]) if len(path)>1)
        self.dimensions=tuple(marks.count(i) for i in self.marked)
        self.ids={}
        self.polynomials=[]
        self.coefficients=defaultdict(lambda:defaultdict(fmpq))
        self.nonzero=[]
        self.work=0
        self.pairs=0

    def exterior(self,poly):
        terms=list(poly.terms())
        assert all(c.denominator==1 for _,c in terms)
        content=gcd(*(int(c) for m,c in terms))
        if terms[0][1]<0:
            content=-content
        key=tuple((m,int(c)//content) for m,c in terms)
        if key not in self.ids:
            self.ids[key]=len(self.polynomials)
            self.polynomials.append(poly/fmpq(content))
        return self.ids[key],content

    def split(self,poly):
        rows=[]
        stars=(self.star,) if self.star else ()
        for deg,local,outside,scale,poles in split_paths(poly,self.marked,self.data[5],stars):
            i,num=self.exterior(outside)
            rows.append((deg,local,i,fmpq(str(scale))*num,sum(poles)))
        return rows

    def simplex_coefficients(self,hard,weighted,weight):
        if self.star or not os.environ.get('MIXED_SIMPLEX_GRAM_LIB'):return False
        from simplex_gram_moments import build
        gram=build(hard,weighted,self.dimensions)
        if gram is None:return False
        denominator,rows=gram
        hd=lcm(*(int(r[3].denominator) for r in hard))
        wd=lcm(*(int(r[3].denominator) for r in weighted))
        h=[int(r[3]*hd) for r in hard]
        w=[int(r[3]*wd) for r in weighted]
        scale=fmpq(str(weight))/(2*denominator*hd**2*wd)
        @cache
        def degrees(a,b):
            return [tuple(k+l+m+u for k,l,m,u in zip(a,b,r[0],self.dimensions)) for r in weighted]
        @cache
        def restore(degree):
            return scale*prod(map(factorial,degree))
        # Collect one left hard row at a time before reinstating the simplex
        # factorials. Products have integral content after the two scale lifts.
        for i,j,values in rows:
            a,b=hard[i],hard[j]
            prefactor=h[i]*h[j]*(1+(i!=j))
            self.pairs+=1
            for k,(degree,moment) in enumerate(zip(degrees(a[0],b[0]),values)):
                if moment:
                    ids=tuple(sorted((a[2],b[2],weighted[k][2])))
                    self.coefficients[degree][ids,0]+=restore(degree)*(prefactor*w[k]*moment)
                    self.work+=1
        return True

    def add(self,color,weight,checkpoints,collect_root=None):
        probe=PrimitiveChamber(self.core_index,self.marks,color)
        groups=[]
        alternate=copy(probe)
        _,alternate.sheared,alternate.flow=sparse_flow(probe.graph,color,probe.sheared)
        if self.adapted:
            from mixed_current_basis import adapted_roots
            roots=adapted_roots(alternate,self.data)
        else:
            roots=((d,scale*self.soft[rays]) for d,rays,scale in square_roots(alternate,self.soft))
        for d,weight_poly in roots:
            for row in groups:
                if d==row[0]:
                    row[1]+=weight_poly
                    break
            else:
                groups.append([d,weight_poly])
        if not groups:
            return
        self.nonzero.append(color)
        assert probe.weight==weight
        for original,old,scale,_,check in checkpoints:
            check[0]=(check[0]+probe.numerator(original,check=True)*old*scale)%PRIME
        t=self.data[1].gens()
        voltage=sum(k*prod(t[i] for i in ex) for ex,k in probe.voltage_terms)
        for root,soft in groups:
            hard=self.split(root)
            weighted=self.split(soft*voltage)
            if self.simplex_coefficients(hard,weighted,weight):
                if collect_root is not None:
                    collect_root(self)
                    self.coefficients.clear()
                    self.ids.clear()
                    self.polynomials.clear()
                continue
            moments=PositionMoments([r[1] for r in weighted])
            factors=[moments.factor(r[1]) for r in hard]
            for i,j in combinations_with_replacement(range(len(hard)),2):
                d,x,h,c,s=hard[i]
                e,y,k,b,u=hard[j]
                degree=tuple(x+y+z for x,y,z in zip(d,e,self.dimensions))
                factor=c*b*fmpq(str(weight*F(1+(i!=j),2)))
                self.pairs+=1
                for (v,_,l,a,w),moment in zip(weighted,moments.product(x,y,factors[i],factors[j]),strict=True):
                    if not moment:
                        continue
                    final=tuple(x+y for x,y in zip(degree,v))
                    pole=11-s-u-w if self.star else 0
                    assert pole>=0
                    key=(tuple(sorted((h,k,l))),pole)
                    self.coefficients[final][key]+=factor*a*moment
                    self.work+=1
            if collect_root is not None:
                collect_root(self)
                self.coefficients.clear()
                self.ids.clear()
                self.polynomials.clear()

    def projected(self):
        star_sum=sum(T[i] for i in self.star)
        jac=prod(T[i] for i in self.star)
        @cache
        def power(n):
            return star_sum**n
        exact=None
        mode=os.environ.get('MIXED_KRONECKER_EXTERIOR')
        if not self.star and self.coefficients and mode in ('1','integer'):
            from mixed_exterior_product import ExteriorProduct
            exact=ExteriorProduct(self.polynomials,self.marked,CTX)
        for degree,terms in self.coefficients.items():
            if exact is not None:
                exterior=exact.contract(terms,mode)
                if exterior is not None:
                    scale=fmpq(1,prod(map(factorial,degree)))
                    yield exterior*scale*prod(T[i]**d for i,d in zip(self.marked,degree))
                    continue
            poles=defaultdict(list)
            for (factors,pole),coefficient in terms.items():
                if coefficient:poles[pole].append((factors,coefficient))
            poly=CTX.constant(0)
            for pole,products in poles.items():
                # Factor common exterior polynomials before expanding the star power.
                frequency=Counter(i for factors,_ in products for i in set(factors))
                def order(factors):
                    return tuple(sorted(factors,key=lambda i:(-frequency[i],i)))
                arranged=sorted((order(factors),coefficient) for factors,coefficient in products)
                block=CTX.constant(0)
                for first,outer in groupby(arranged,key=lambda r:r[0][0]):
                    inner=CTX.constant(0)
                    for second,leaves in groupby(outer,key=lambda r:r[0][1]):
                        last=sum(c*self.polynomials[f[2]] for f,c in leaves)
                        inner+=self.polynomials[second]*last
                    block+=self.polynomials[first]*inner
                poly+=block*power(pole)
            if self.star:
                poly,remainder=divmod(poly,star_sum)
                assert not remainder
            scale=fmpq(1,prod(map(factorial,degree)))
            yield poly*jac*scale*prod(T[i]**d for i,d in zip(self.marked,degree))


def direct_checks(marks,star):
    checks=[]
    for metric in (METRIC,(11,5,2,7,9,6,8,3,4)):
        probe=PrimitiveChamber(3,marks)
        new=list(map(F,metric))
        delta=sum(new[i]*new[j] for i,j in combinations(star,2))
        for i in star:
            new[i]=delta/new[i]
        jac=prod(new[i] for i in star)/sum(new[i] for i in star)**3
        tree=kirchhoff(reduction(marks,star)[0],new)
        checks.append((metric,pow(probe.kirchhoff(metric),-13,PRIME),residue(jac*tree**13),
                       tuple(fmpq(str(x)) for x in new),[0,F()]))
    return checks


def projected_placement(original,original_star,words=None,target=(0,1,2)):
    """Combine any set of hard orders using the same empty prism vertex."""
    probe=PrimitiveChamber(3,original)
    move=next(g for g in probe.actions if tuple(sorted(g[i][0] for i in original_star))==target)
    m=tuple(sorted(move[i][0] for i in original))
    s=target
    batch=ColorMoments(m,s)
    checks=direct_checks(m,s)
    for color,n in chamber_rows(original):
        if words is not None and ''.join(map(str,color)) not in words:
            continue
        config=tuple(tuple(color[i] for i in row) for row in probe.groups)
        cc=sum(image(config,move),())
        batch.add(cc,F(n,12),checks)
    def slices():
        for poly in batch.projected():
            for _,_,_,metric,totals in checks:
                totals[1]+=F(str(poly(*metric)))
            yield poly
        for _,_,_,_,totals in checks:
            assert residue(totals[1])==totals[0]
    bank,count=wheel_bank(slices(),m,s)
    bins=sum(sum(bool(c) for c in terms.values()) for terms in batch.coefficients.values())
    stats=[len(batch.nonzero),batch.pairs,batch.work,bins,len(batch.polynomials),count,len(bank)]
    return bank,stats

def geometry():
    core=PrimitiveChamber(3).core
    triangles=((0,1,3),(6,7,8))
    actions=PrimitiveChamber(3).actions
    population=defaultdict(set)
    for rung in (2,4,5):
        for first,other in permutations(triangles):
            for double,single in product(first,other):
                if set(core[rung]).isdisjoint(core[double]):
                    continue
                marks=tuple(sorted((rung,)*3+(double,)*2+(single,)))
                canonical=min(tuple(sorted(g[i][0] for i in marks)) for g in actions)
                population[canonical].add(marks)
    assert population.keys()=={(0,0,2,2,2,i) for i in (6,7,8)}
    assert all(len(shapes)==12 for shapes in population.values())
    return population


def control():
    marks=(0,0,2,2,2,6)
    words=('001122','001212','002112','010122','010212')
    joint,counts=projected_placement(marks,ORIGINAL_STAR,words)
    separate=defaultdict(F)
    for word in words:
        color=tuple(map(int,word))
        bank,_=wheel_bank(transformed(color,True,marks,ORIGINAL_STAR),marks,ORIGINAL_STAR)
        for k,v in bank.items():
            separate[k]+=v
    assert joint=={k:v for k,v in separate.items() if v}
    assert counts==[5,15,525,126,32,26278,13332]
    print('Independent color projections agree on',len(joint),'four-wheel powers.')


def exact_columns(record):
    cols=tuple(map(F,record['columns']))
    lo,hi=bounds(cols)
    assert record['bound']<10**20*lo<=10**20*hi<record['bound']+1
    return cols


def replay(records,selected,rebuild):
    totals={}
    banks={}
    for marks in sorted(geometry()):
        edge=marks[-1]
        rec=records['placements'][str(edge)]
        if selected is not None and selected!=edge:
            continue
        totals[edge]=exact_columns(rec)
        banks[edge]=defaultdict(fmpq)
        rows={''.join(map(str,c)):n for c,n in chamber_rows(marks)}
        assert len(rows)==90 and set(rows.values())=={96}
        coverage=set()
        assert sum(stats[0] for _,_,stats in rec['blocks'])==54
        for prefix,center,stats in rec['blocks']:
            words={word for word in rows if word.startswith(prefix)}
            assert words and coverage.isdisjoint(words)
            coverage|=words
            if rebuild:
                star=tuple(i for i,e in enumerate(PrimitiveChamber(3,marks).core) if center in e)
                # Reap the multivariate determinant workspace after each block.
                with get_context('spawn').Pool(1) as worker:
                    part,actual=worker.apply(projected_placement,(marks,ORIGINAL_STAR,words,star))
                assert actual==stats
                for key,v in part.items():
                    banks[edge][key]+=fmpq(str(v))
        assert coverage==rows.keys()
        if rebuild:
            assert sum(bool(v) for v in banks[edge].values())==rec['wheel_orbits']
        print(marks,'54 nonzero colors;',rec['wheel_orbits'],'wheel powers;',*totals[edge],flush=True)
    if rebuild:
        with TemporaryDirectory() as tmp:
            keys=set.union(*(set(b) for b in banks.values()))
            values=wheel_values(keys,Path(tmp))
            for edge,b in banks.items():
                assert capture(b,values)==tuple(fmpq(str(v)) for v in totals[edge])
    if 7 in totals:
        two_star=json.loads((HERE/'mixed_prism_empty_stars.json').read_text())
        assert totals[7]==(*map(F,two_star['231']['sum']),F())
    if selected is None:
        for name,edges in (('sum',(6,7,8)),('new',(6,8))):
            assert tuple(map(sum,zip(*(totals[i] for i in edges))))==exact_columns(records[name])
        weight_five=exact_columns(records['new'])[2]
        for table in ('mixed_prism_one_star.json','mixed_prism_opposite_face.json'):
            weight_five+=F(json.loads((HERE/table).read_text())['sum'][2])
        assert F(3,16)*weight_five==-F(8818959362641,129060195264000)
        print('36 placements, 162 nonzero color orbits, 15552 labeled orders:',*exact_columns(records['sum']))


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true',help='Compare joint and independent color projections')
    parser.add_argument('--rebuild',action='store_true',help='Recalculate the signed moment tensors and angular periods')
    parser.add_argument('--edge',type=int,choices=(6,7,8),help='Select the singleton path')
    opts=parser.parse_args()
    records=json.loads(Path(__file__).with_suffix('.json').read_text())
    replay(records,opts.edge,opts.rebuild)
    if opts.check:
        control()
