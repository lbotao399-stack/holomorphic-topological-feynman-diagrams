"""Contract available perfect-matching color polynomials against the exact dimension-26 integrals."""
from argparse import ArgumentParser
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
from time import sleep
import gzip
import json

from mixed_bare_symbolic import classified_orbits
from mixed_bipartite_star_bank import measure
from mixed_bipartite_three_paths import bounds


RESULT=Path(__file__).with_suffix('.json')
WORDS=('005588','000558')


def digest(source):
    h=sha256()
    with gzip.open(source,'rb') as data:
        for block in iter(lambda:data.read(1048576),b''):h.update(block)
    return h.hexdigest()


def columns(row):
    return tuple(map(Fraction,row['period']))


def package(v):
    v=tuple(map(Fraction,v))
    return {'period':list(map(str,v)),'bounds':bounds(v),'ward_zeta5':str(Fraction(3,16)*v[2])}


def run(scratch,wait):
    bank=scratch/'matching13'
    records=json.loads(RESULT.read_text()) if RESULT.exists() else {}
    rows=classified_orbits(4)[0]
    cases={w:[(colour,count) for m,colour,count,_ in rows if m==tuple(map(int,w))] for w in WORDS}
    def save():
        pending=RESULT.with_suffix('.partial')
        pending.write_text(json.dumps(records,sort_keys=True,indent=2)+'\n')
        pending.replace(RESULT)
    checked=set()
    while True:
        for word in WORDS:
            entry=records.setdefault(word,{'colors':{}})
            for color,count in cases[word]:
                tag=''.join(map(str,color))
                source=scratch/f'{word}-{tag}.gz'
                if not source.exists() or (word,tag) in checked:continue
                signature=digest(source)
                if tag not in entry['colors']:
                    v=measure(bank,13,scratch,source,10000)
                    entry['colors'][tag]=dict(package(v),nonzero_hard_orders=count,weights_sha256=signature)
                    save()
                    print(word,tag,*v,flush=True)
                actual=entry['colors'][tag]
                assert actual==dict(package(columns(actual)),nonzero_hard_orders=count,weights_sha256=signature)
                checked.add((word,tag))
            combined=scratch/f'{word}-direct.gz'
            if len(entry['colors'])==len(cases[word]) and combined.exists() and 'total' not in entry:
                v=tuple(sum(c) for c in zip(*(columns(r) for r in entry['colors'].values()),strict=True))
                assert tuple(map(str,v))==measure(bank,13,scratch,combined,10000)
                entry['total']=dict(package(v),nonzero_hard_orders=sum(n for _,n in cases[word]),
                                    weights_sha256=digest(combined))
                save()
                print(word,'complete',entry['total'],flush=True)
        if not wait or all('total' in records[w] for w in WORDS):break
        sleep(3)


def low_occupancy():
    from collections import Counter
    from mixed_bipartite_three_paths import RESULT as THREE_PATHS
    from mixed_bipartite_triples import PERIOD

    matching=json.loads(RESULT.read_text())
    empty=json.loads(THREE_PATHS.read_text())
    by_word={'000555':PERIOD}
    by_word.update({word:columns(row) for word,row in empty.items()})
    by_word.update({word:columns(matching[word]['total']) for word in WORDS})
    prime=1009
    rows=[[v.numerator*pow(v.denominator,-1,prime)%prime for v in by_word[word]]
          for word in ('000558','005588','000113')]
    (a,b,c),(d,e,f),(g,h,i)=rows
    assert (a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g))%prime==474
    survivors=[(m,n) for m,c,n,_ in classified_orbits(4)[0] if len(set(m))<=3]
    counts=Counter(''.join(map(str,m)) for m,n in survivors)
    assert set(by_word)==set(counts)
    records={}
    for label,words in (('three_paths',[w for w in by_word if len(set(w))==3]),
                        ('at_most_three_paths',list(by_word))):
        values=tuple(sum(i) for i in zip(*(by_word[w] for w in words),strict=True))
        records[label]=dict(package(values),placements=sorted(words),
                            nonzero_colors=sum(counts[w] for w in words),
                            nonzero_hard_orders=sum(n for m,n in survivors if ''.join(map(str,m)) in words))
    output=RESULT.with_name('mixed_bipartite_low_occupancy.json')
    pending=output.with_suffix('.partial')
    pending.write_text(json.dumps(records,sort_keys=True,indent=2)+'\n')
    pending.replace(output)
    return records


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('scratch',type=Path)
    parser.add_argument('--wait',action='store_true',help='consume newly completed color and whole-placement weights')
    args=parser.parse_args()
    run(args.scratch,args.wait)
    if all('total' in row for row in json.loads(RESULT.read_text()).values()):low_occupancy()
