"""Evaluate invariant quadratic-soft potentials by S_8 stabilizer Gram matrices."""
from itertools import combinations,combinations_with_replacement
from pathlib import Path
from random import Random
from subprocess import run
from tempfile import TemporaryDirectory
from genus_two_probes import det
from quartic_unisolvence import pivot_rows


SUBSETS=list(combinations(range(10),5))


def candidates(seed,k):
    rng=Random(seed)
    triples=[]
    for i in range(k):
        size=rng.choices((1,2,3),weights=(1,3,8))[0]
        w=[rng.randrange(252) for _ in range(size)]
        w+=rng.choices(w,k=3-size)
        triples.append(tuple(sorted(w)))
    words=[]
    for a,b in combinations_with_replacement(range(k),2):
        w=triples[a]+triples[b]
        flags=[(64 if i==0 else 128 if i==1 else 0)+sum(1<<j for j,s in enumerate(w) if i in SUBSETS[s])
               for i in range(10)]
        if any((h&63).bit_count()+int(i==1)<2 for i,h in enumerate(flags)):
            continue
        if any(flags.count(h)>1 and (h&63).bit_count()%2 for h in flags):
            continue
        words.append((a,b))
    return triples,words


def fixtures(count):
    rng=Random(78965)
    out=[]
    for i in range(count):
        rows=[tuple(int(j==h) for j in range(5)) for h in range(5)]
        rows += [tuple(rng.randrange(-4,5) for h in range(5)) for _ in range(4)]
        rows.append(tuple(-sum(x) for x in zip(*rows)))
        out.append(rows)
    return out


def physical(p,a):
    coeff=[]; gauge=[]
    for i in range(9):
        h=next(j for j,x in enumerate(p[i]) if x)
        for j in range(5):
            if j!=h:
                coeff.append(p[i][h]*(a[5*i+j]-a[45+j])-p[i][j]*(a[5*i+h]-a[45+h]))
                gauge.append(p[i][h]*p[9][j]-p[i][j]*p[9][h])
    h=next(i for i,x in enumerate(gauge) if x%1009)
    return [(gauge[h]*x-coeff[h]*gauge[j])%1009 for j,x in enumerate(coeff) if j!=h]


def evaluate(exe,points,triples,words):
    data=[(len(points),len(triples),len(words)),*triples,*words]
    for p in points:
        data.extend(p)
        data.append([det([p[i] for i in t])%1009 for t in SUBSETS])
    stdin='\n'.join(' '.join(map(str,r)) for r in data)
    result=run([exe],input=stdin,text=True,capture_output=True,check=True)
    lines=iter(result.stdout.splitlines())
    out=[]
    for p in points:
        cols=[physical(p,list(map(int,next(lines).split()))) for w in words]
        out.extend(zip(*cols))
    return out


def basis_data():
    import json
    meta=json.loads(Path(__file__).with_name('soft_target_basis.json').read_text())
    triples=sorted({tuple(t) for w in meta['words'] for t in (w[:3],w[3:])})
    indices={t:i for i,t in enumerate(triples)}
    words=[tuple(sorted((indices[tuple(w[:3])],indices[tuple(w[3:])]))) for w in meta['words']]
    return triples,words,meta


def check(exe):
    from flint import nmod_mat
    triples,words,meta=basis_data()
    rows=evaluate(exe,fixtures(15),triples,words)
    selected=[row for i,row in enumerate(rows) if i not in meta['omitted_rows']]
    assert len(selected)==len(words)==518
    assert int(nmod_mat(selected,1009).det())==meta['determinant']==356
    print('518 soft potential classes tested on 15 fixtures; determinant 356 modulo 1009.')


if __name__=='__main__':
    with TemporaryDirectory() as tmp:
        exe=Path(tmp)/'soft-target-orbits'
        run(['c++','-O3','-std=c++17','-DACCELERATE_NEW_LAPACK','-framework','Accelerate',
             str(Path(__file__).with_name('soft_target_orbits.cpp')),'-o',str(exe)],check=True)
        check(exe)
