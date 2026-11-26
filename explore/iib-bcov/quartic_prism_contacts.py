"""Integrated prism dimension shift including all three rung contacts."""
from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
from tempfile import TemporaryDirectory
import json

from quartic_prism_shift import symbols
from quartic_projective_ibp import faces, fold, PERMS, PRISM, reduction
from quartic_terminal_periods import cut_periods, WHEEL_PERMS


def controls(folder,executable):
    u,a,b,d,_,_=symbols()
    x=u.context().gens()
    rungs=(0,6,7)
    edges=range(9)
    cotrees=[tuple(map(int,t)) for t,c in u.terms()]
    subsets=[(ids,min(sum(t[i] for i in ids) for t in cotrees)) for n in range(1,9)
             for ids in combinations(edges,n)]
    def check(p,q,removed=None):
        # Restricting U to t_e=0 selects cotrees not containing that edge.
        for ids,rank in subsets:
            if removed in ids:continue
            if removed is not None:
                if len(ids)==8:continue
                rank=min(sum(t[i] for i in ids) for t in cotrees if not t[removed])
            assert sum(p[i] for i in ids)>q*rank
    def prism(p,q):
        check(p,q)
        return fold(p,PERMS)
    def wheel(p,e):
        check(p,3,e)
        _,contracted=faces(PRISM,e)
        moves,name,keep=reduction(contracted)
        assert name=='wheel' and not moves
        shortened=tuple(p[i] for i in edges if i!=e)
        return fold(tuple(shortened[j] for j in keep),WHEEL_PERMS)
    measurements=[]
    for i,j in product((1,2,5),(3,4,8)):
        m=tuple(int(e in (*rungs,i,j)) for e in edges)
        h=u.context().from_dict({m:1})
        left=Counter()
        for ex,c in (h*a*b).terms():
            left[prism(tuple(int(k)+1 for k in ex),4)]+=int(c)
        interior=Counter()
        contact=Counter()
        for e,f in combinations(rungs,2):
            p=tuple(k+1-int(v in (e,f)) for v,k in enumerate(m))
            interior[prism(p,3)]+=2
        for e in rungs:
            p=tuple(k+1 for k in m)
            contact[wheel(p,e)]-=1
        measurements.append((left,interior,contact))
    path=Path(__file__).with_name('quartic_terminal_periods.json')
    terminal=json.loads(path.read_text())
    sources=sorted(set().union(*(set(left) for left,_,_ in measurements)))
    upper=cut_periods('prism',sources,4,folder,executable)
    def bank(rows):
        return {tuple(r['powers']):tuple(map(F,r['1,zeta3,zeta5'])) for r in rows}
    tables=[bank(upper),bank(terminal['prism']),bank(terminal['wheel'])]
    results=[]
    for i,(left,bulk,faces_) in enumerate(measurements):
        def period(row,table):
            return tuple(sum(v*table[p][j] for p,v in row.items()) for j in range(3))
        high,interior,contact=(period(row,table) for row,table in zip((left,bulk,faces_),tables))
        assert tuple(30*v for v in high)==tuple(b+c for b,c in zip(interior,contact))
        assert any(contact)
        results.append((high,interior,contact))
    print('q=4 to q=3 including all rung contacts:',len(measurements),'measurements,',len(sources),
          'upper moment orbits,',len(set(results)),'distinct integrated triples',flush=True)
    for result in dict.fromkeys(results):print(result,flush=True)
    return results


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('form',type=Path)
    parser.add_argument('forcer',type=Path)
    args=parser.parse_args()
    with TemporaryDirectory(prefix='bcov-rung-contact-') as temporary:
        folder=Path(temporary)
        for item in ('forcer','forcer.h'):
            (folder/item).symlink_to((args.forcer/item).resolve(strict=True))
        controls(folder,args.form.resolve(strict=True))
