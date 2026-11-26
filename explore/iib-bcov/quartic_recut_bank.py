"""Enumerate the planar and nine-cut residue functionals on the physical K3,3 bank."""
from collections import Counter,defaultdict
from functools import cache
from gzip import open as open_gzip
from itertools import islice
from pathlib import Path
import numpy as np
from quartic_nonplanar_cut import AUTOMORPHISMS,endpoint_gram,homogeneous_endpoint
from quartic_planar_vertex import recut_permutations


@cache
def endpoint_tables():
    groups=defaultdict(list)
    for g in AUTOMORPHISMS:
        groups[g[0],frozenset((g[6],g[7]))].append(g)
    keys=sorted(groups,key=lambda key:(key[0],sorted(key[1])))
    assert len(keys)==18 and all(len(groups[key])==4 for key in keys)
    weights=((0,1,0,1,0,1,1,1),(1,0,1,0,1,1,1,0),
             (0,1,1,0,1,2,0,1),(1,0,0,1,1,0,2,1))
    for w in weights:
        assert all(sum(a*b for a,b in zip(w,e[:8]))+e[8]>=2 for e,c in endpoint_gram(0).terms())
    lower=[]
    for j in range(9):
        g=next(g for g in AUTOMORPHISMS if g[0]==j)
        vertices=[(g[r[0]],frozenset((g[r[6]],g[r[7]]))) for r in recut_permutations()]
        lower.append([keys.index(v) for v in vertices])
    assert list(sorted(Counter(i for row in lower for i in row).values()))==[2]*18
    return [groups[k] for k in keys],np.array(lower),np.array(weights,dtype=np.int16)


def physical_rows(independent=False):
    groups,lower,tests=endpoint_tables()
    code=16**np.arange(7,-1,-1,dtype=np.int64)
    scores=[np.array([[0 if k==g[0] else int(code[g.index(k)-1]) for k in range(9)]
                      for g in group],dtype=np.int64).T for group in groups]
    def scaled(word):
        num,sep,den=word.partition('/')
        denominator=int(den or 1)
        assert 1925%denominator==0
        return int(num)*(1925//denominator)
    with open_gzip(Path(__file__).with_name('quartic_k33.weights.gz'),'rt') as stream:
        while lines:=list(islice(stream,2048)):
            words=[row.split() for row in lines]
            alpha=np.array([r[:9] for r in words],dtype=np.int16)-1
            weights=np.array([[scaled(v) for v in row[9:]] for row in words],dtype=np.int64)
            alive=np.empty((len(alpha),18),np.bool_)
            keys=np.empty((len(alpha),18),np.int64)
            for i,group in enumerate(groups):
                uncut=alpha[:,group[0][1:]]
                alive[:,i]=(uncut@tests.T>=20).all(axis=1)
                keys[:,i]=(alpha@scores[i]).min(axis=1)
            zero=~alive.reshape(len(alpha),9,2).any(axis=2)
            allowed=np.where(zero.any(axis=1)[:,None],zero,True)
            # A highest-power propagator and then the least ordinary uncut
            # representative fixes both previous maximal-cut populations.
            levels=np.where(allowed,alpha,-1)
            allowed &= levels==levels.max(axis=1)[:,None]
            ordinary=keys.reshape(len(alpha),9,2).min(axis=2)
            choice=np.where(allowed,ordinary,2**32).argmin(axis=1)
            if independent:
                # Re-evaluate both endpoints of each disjoint cut. Reversing the
                # canonical endpoint representative independently checks label covariance.
                keys=np.stack([(alpha@matrix).max(axis=1) for matrix in scores],axis=1)
                lower_used=np.array([[endpoint for edge in sorted(set(v//2 for v in row))
                                      for endpoint in (2*edge,2*edge+1)] for row in lower])
            else:lower_used=lower
            for i in range(len(alpha)):
                yield (keys[i,alive[i]],keys[i,lower_used[choice[i]]][alive[i,lower_used[choice[i]]]],
                       weights[i],bool(zero[i].any()))


def summarize():
    zeros,n,top,planar=0,0,Counter(),Counter()
    for all_keys,selected,weights,zero in physical_rows():
        zeros+=int(zero)
        n+=1
        def acc(keys,hist):
            for k in keys:
                degree=sum((int(k)>>(4*j))&15 for j in range(8))-20
                hist[degree]+=1
        acc(all_keys,top)
        acc(selected,planar)
    assert (n,zeros)==(465564,395810)
    print('degree counts for the nine-cut sum',sorted(top.items()),flush=True)
    print('degree counts for the planar double pole',sorted(planar.items()),flush=True)


def native_input(output,independent=False):
    from struct import Struct
    row=Struct('<IqqB')
    n=zero_count=0
    with output.open('wb') as stream:
        for all_keys,selected,weights,zero in physical_rows(independent):
            chosen=Counter(map(int,selected))
            for key,count in Counter(map(int,all_keys)).items():
                assert count<=255 and chosen[key]<=15 and count<=15
                stream.write(row.pack(key,*weights,count+16*chosen[key]))
            n+=1
            zero_count+=int(zero)
    assert (n,zero_count)==(465564,395810)


def compare_native(binary,ops):
    from fractions import Fraction as F
    from math import factorial,prod
    from struct import Struct
    from tempfile import TemporaryDirectory
    import subprocess
    row=(9,6,5,6,4,5,4,5,4)
    norm=F(prod(factorial(v-1) for v in row),factorial(11))
    expected_top=F(1693,21339318)/(20*norm)
    groups,lower,tests=endpoint_tables()
    with TemporaryDirectory(prefix='bcov-recut-control-') as folder:
        folder=Path(folder)
        (folder/'parts').mkdir()
        with (folder/'weights').open('wb') as stream:
            for i,group in enumerate(groups):
                ex=min(tuple(row[j]-1 for j in g[1:]) for g in group)
                key=sum(int(v)<<(28-4*j) for j,v in enumerate(ex))
                stream.write(Struct('<IqqB').pack(key,1925,-3850,1+16*int(i in lower[0])))
        output=subprocess.run([binary,ops,folder/'weights',folder/'parts'],check=True,capture_output=True,text=True)
        values=tuple(map(F,output.stdout.split()))
        assert values==(expected_top,-2*expected_top,F(1693,37072792505548800),-F(1693,18536396252774400))
        print(output.stderr,flush=True)
        print('Physical nine-cut sum times its projective weight and planar double pole',
              values[0]*20*norm,values[2],flush=True)
    return values


def terminal_controls():
    import json
    from fractions import Fraction as F
    from math import factorial,prod
    from quartic_nonplanar_cut import projection
    from quartic_projective_ibp import K33,faces,reduction
    rows=json.loads(Path(__file__).with_name('quartic_terminal_periods.json').read_text())
    nonplanar=[next(g for g in AUTOMORPHISMS if g[0]==i) for i in range(9)]
    moves,name,keep=reduction(faces(K33,5)[1])
    assert not moves and name=='wheel'
    for graph in ('K33','wheel'):
        for row in rows[graph]:
            source=row['powers']
            if graph=='wheel':
                mapped=[0]*8
                for i,j in enumerate(keep):mapped[j]=source[i]
                full=mapped[:5]+[0]+mapped[5:]
            else:full=source
            weight=F(prod(factorial(a-1) for a in source),2)
            total=sum(projection(tuple(full[i] for i in g[1:]),1) for g in nonplanar)
            assert 20*weight*total==F(row['1,zeta3,zeta5'][2])
    print('All twelve terminal bipartite and wheel periods pass the nine-cut weight-five test.',flush=True)


def finite_controls(form,forcer,q):
    from fractions import Fraction as F
    from itertools import combinations,product
    from math import factorial,prod
    from tempfile import TemporaryDirectory
    from quartic_nonplanar_cut import projection
    from quartic_planar_periods import forest_data
    from quartic_projective_ibp import K33,fold
    from quartic_terminal_periods import cut_periods
    terminal_controls()
    tree=forest_data(tuple(K33),6)[1]
    cotrees=[tuple(map(int,e)) for e,c in tree.terms()]
    bounds=[(I,min(sum(t[i] for i in I) for t in cotrees)) for k in range(1,9)
            for I in combinations(range(9),k)]
    rows=[p for p in product(range(1,q),repeat=9) if sum(p)==4*q and p==fold(p,AUTOMORPHISMS)
          and all(sum(p[i] for i in I)>q*h for I,h in bounds)]
    assert q!=4 or len(rows)==35
    with TemporaryDirectory(prefix='bcov-nine-cut-controls-') as folder:
        folder=Path(folder)
        for name in ('forcer','forcer.h'):(folder/name).symlink_to((forcer/name).resolve(strict=True))
        reference=cut_periods('K33',rows,q,folder,form.resolve(strict=True))
    group=[next(g for g in AUTOMORPHISMS if g[0]==i) for i in range(9)]
    for row,p in zip(reference,rows):
        assert tuple(row['powers'])==p
        weight=F(prod(factorial(a-1) for a in p),factorial(q-1))
        total=sum(projection(tuple(p[i] for i in g[1:]),q-2) for g in group)
        actual=20*weight*total
        expected=F(row['1,zeta3,zeta5'][2])
        assert actual==expected,(p,actual,expected)
    print('Nine top cuts reproduce the full zeta(5) coefficient on all',len(rows),'finite K3,3 orbits at q =',q,flush=True)
    return reference


def write_operators(path):
    with path.open('w') as stream:
        for terms in homogeneous_endpoint(0)[1]:
            print(len(terms),file=stream)
            for ex,c in terms:print(*ex,c,file=stream)


def certificate(binary):
    from fractions import Fraction as F
    from subprocess import run
    from tempfile import TemporaryDirectory
    expected=(F(249372365272156031,16132524408000000),F(2209968587039007449,27877002177024000000),
              -F(169378632737740819356497,149403308542488000000),-F(746978483911572523339,241054077648384000000))
    with TemporaryDirectory(prefix='bcov-nine-cut-bank-') as folder:
        folder=Path(folder)
        ops=folder/'operators'
        write_operators(ops)
        compare_native(binary,ops)
        weights=folder/'weights'
        for independent in (False,True):
            native_input(weights,independent)
            with TemporaryDirectory(dir=folder) as pages:
                result=run([binary,ops,weights,pages],capture_output=True,text=True,check=True)
            values=tuple(map(F,result.stdout.split()))
            if independent:values=values[:2]+tuple(v/2 for v in values[2:])
            assert values==expected
            print(result.stderr,flush=True)
    print('Nine-cut sums and planar double poles on the complete K3,3 bank:',*expected,flush=True)
    return expected


if __name__=='__main__':
    from argparse import ArgumentParser
    p=ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path)
    p.add_argument('--operators',type=Path)
    p.add_argument('--independent',action='store_true')
    p.add_argument('--compare',nargs=2,type=Path)
    p.add_argument('--finite',nargs=2,type=Path)
    p.add_argument('--verify',type=Path)
    args=p.parse_args()
    if args.verify:certificate(args.verify.resolve(strict=True))
    elif args.finite:finite_controls(*args.finite,4)
    elif args.compare:compare_native(*args.compare)
    elif args.input:native_input(args.input,args.independent)
    elif args.operators:write_operators(args.operators)
    else:summarize()
