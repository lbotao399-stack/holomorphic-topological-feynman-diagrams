"""Finish the mixed bipartite one-star chamber by convergent six-point identities."""
from argparse import ArgumentParser
from collections import defaultdict
from fractions import Fraction as F
from math import factorial,prod
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory
import gzip
import json

from mixed_bipartite_fixed_spin import (
    HERE,compile_program,evaluate,reflected,verify_chart,zeta_interval,
)


def reduced(seeds,q,ibp,wheel,tmp,wheel_cache=None,prime_workers=1,wheel_reference=None,prism_keys=None,
            return_leaves=False):
    output=tmp/'reduction'
    rows=''.join(' '.join(map(str,seed))+'\n' for seed in seeds)
    command=[ibp,str(q),output]
    if prism_keys is not None:command.append(prism_keys)
    run(command,input=rows,text=True,check=True)
    original={}
    pivots={}
    leaves={}
    wheels=defaultdict(set)
    for line in output.read_text().splitlines():
        tag,*row=line.split()
        if tag=='I':
            i,j,z=row
            original[int(i)]=(int(j),F(z))
        elif tag=='P':
            i,j,z=row
            pivots.setdefault(int(i),[]).append((int(j),F(z)))
        else:
            assert tag=='T'
            i,d,k,z=row
            d,k=int(d),int(k)
            leaves[int(i)]=(d,k,F(z))
            wheels[d].add(k)
    assert len(original)==len(seeds)
    if return_leaves:
        expansions={i:{(d,k):scale} for i,(d,k,scale) in leaves.items()}
        for i,terms in reversed(pivots.items()):
            combined=defaultdict(F)
            for j,scale in terms:
                for key,weight in expansions[j].items():combined[key]+=scale*weight
            expansions[i]={k:c for k,c in combined.items() if c}
        return {s:{k:z*c for k,c in expansions[j].items()} for s,(j,z) in original.items()}
    table={}
    for d,keys in sorted(wheels.items()):
        known={}
        if d==-1:
            for db in (wheel_cache,wheel_reference):
                if db is None:continue
                wanted=list(keys-known.keys())
                for start in range(0,len(wanted),500):
                    batch=wanted[start:start+500]
                    marks=','.join('?' for _ in batch)
                    for k,r,w3,w5 in db.execute(f'SELECT key,r,w3,w5 FROM prism WHERE key IN ({marks})',batch):
                        known[k]=tuple(map(F,(r,w3,w5)))
            if known.keys()!=keys:raise ArithmeticError('Incomplete angular prism reference')
            table[d]=known
            continue
        for db in (wheel_cache,wheel_reference):
            if db is None:continue
            assert q==13
            wanted=list(keys-known.keys())
            for start in range(0,len(wanted),500):
                batch=wanted[start:start+500]
                marks=','.join('?' for _ in batch)
                query=f'SELECT key,r,w3,w5 FROM wheels WHERE degree=? AND key IN ({marks})'
                for k,r,w3,w5 in db.execute(query,(d,*batch)):
                    known[k]=tuple(map(F,(r,w3,w5)))
        missing=keys-known.keys()
        for k in missing:verify_chart(k,q,d)
        if missing:
            new=evaluate(q,d,missing,wheel,tmp,prime_workers=prime_workers)
            if wheel_cache is not None:
                wheel_cache.executemany('INSERT OR IGNORE INTO wheels VALUES (?,?,?,?,?)',
                                        ((d,k,*map(str,v)) for k,v in new.items()))
                wheel_cache.commit()
            known.update(new)
        table[d]=known
    values={i:tuple(scale*x for x in table[d][k]) for i,(d,k,scale) in leaves.items()}
    for i,terms in reversed(pivots.items()):
        values[i]=tuple(sum(scale*values[j][c] for j,scale in terms) for c in range(3))
    return {s:tuple(z*x for x in values[j]) for s,(j,z) in original.items()}


def control(ibp,wheel,tmp):
    ref=json.loads((HERE.parents[1]/'nonSP/verify/prism_forcer/periods.json').read_text())
    seeds=[tuple(a+1 for a in ex) for ex in ref['exponents']]
    values=reduced(seeds,6,ibp,wheel,tmp)
    for i,(ex,cols) in enumerate(zip(ref['exponents'],ref['coefficients'],strict=True)):
        assert values[i]==tuple(F(str(v))*factorial(5)/prod(map(factorial,ex)) for v in cols)
    print('Independent q=6 three-column propagators:',len(seeds),flush=True)


KNOWN=(F(1875840319482653063081961833881,8933882473602043013726208000000000),
       F(192571615355536381505743,342912864219655137408000000),
       -F(119708316672311,140256167203152000))


CHAMBER=(F(213249810833125731948743,733743238248244224000000000),
         F(17962743011421743,23108567594112000000),
         -F(859794979,727783056000))


def entire_chamber(ibp,wheel,tmp):
    with gzip.open(HERE/'mixed_bipartite_one_star.components.gz','rt') as stream:
        rows=list(reversed([line.split() for line in stream]))
    seeds=[tuple(map(int,r[:9])) for r in rows]
    assert len(seeds)==len(set(seeds))==37391
    values=reduced(seeds,13,ibp,wheel,tmp)
    total=tuple(sum(F(r[9])*values[i][c] for i,r in enumerate(rows)) for c in range(3))
    assert total==CHAMBER
    print('Direct three-column chamber check on all original seeds:',*total,flush=True)


def physical(complete,ibp,wheel,tmp):
    with gzip.open(HERE/'mixed_bipartite_one_star.components.gz','rt') as stream:
        data={tuple(map(int,r[:9])):r[:11] for line in stream if (r:=line.split())[11]=='0'}
    for kind in ('wheel',1,2,3,4,5):
        for seed in reflected(complete,13,data.values(),kind):del data[seed]
    for d in range(1,7):
        for seed in reflected(complete,13,data.values(),d,formal=True):del data[seed]
    assert len(data)==3840
    seeds=list(data)
    values=reduced(seeds,13,ibp,wheel,tmp)
    add=tuple(sum(F(data[s][9])*values[i][c] for i,s in enumerate(seeds)) for c in range(3))
    total=tuple(a+b for a,b in zip(add,KNOWN,strict=True))
    assert total==CHAMBER
    print('3840 residual components:',*add,flush=True)
    print('complete 37391-component color chamber:',*total,flush=True)
    lo=hi=total[0]
    for v,n in zip(total[1:],(3,5),strict=True):
        a,b=sorted(v*z for z in zeta_interval(n,terms=96))
        lo+=a
        hi+=b
    assert -206524<10**20*lo<=10**20*hi<-206523
    print('strict bounds, units of 10^-20:',(10**20*lo)//1,-((-10**20*hi)//1),flush=True)


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true')
    parser.add_argument('--controls-only',action='store_true')
    parser.add_argument('--entire',action='store_true')
    args=parser.parse_args()
    with TemporaryDirectory() as directory:
        tmp=Path(directory)
        complete,ibp,wheel=tmp/'complete',tmp/'ibp',tmp/'wheel'
        compile_program('quartic_prism_completion.cpp',complete)
        compile_program('quartic_prism_local_ibp.cpp',ibp)
        compile_program('quartic_wheel_bank.cpp',wheel)
        if args.check or args.controls_only:control(ibp,wheel,tmp)
        if args.entire:entire_chamber(ibp,wheel,tmp)
        elif not args.controls_only:physical(complete,ibp,wheel,tmp)
