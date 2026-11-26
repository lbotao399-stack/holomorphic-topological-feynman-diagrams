"""Re-sum the 47 bare cyclic graph records before applying the period functional."""
from collections import defaultdict
from fractions import Fraction
from functools import cache
import json

from cyclic_series import ROOT, SOURCE, engine


def check():
    graphs=json.loads((ROOT/'deepen/o3/higher_orders/four-loop-graphs.json').read_text())
    all_indices=[]
    for core in ('core03-10','core11-38','core39-49'):
        records=json.loads((SOURCE/(core+'.json')).read_text())
        _, edges=records[0]['graphs'],records[0]['edges']
        group=engine(records[0]).cotree_permutations()
        @cache
        def orbit(exponent):
            return min(tuple(exponent[i] for i in permutation) for permutation in group)
        total=defaultdict(lambda:defaultdict(Fraction))
        first,last=records[0]['graphs']
        for i in range(first,last+1):
            rows=json.loads((SOURCE/f'graph{i:02}.json').read_text())
            all_indices.append(i)
            for row in rows:
                assert row['graph']==i and row['edges']==edges
                assert row['automorphism_order']==graphs[i]['aut']
                for exponent,c in row['terms']:
                    full=tuple(a+b for a,b in zip(exponent,row['h']))
                    total[row['N_power']][orbit(full)]+=Fraction(c,row['denominator'])
        total={k:{a:c for a,c in by_exp.items() if c} for k,by_exp in total.items()}
        expected={r['N_power']:{tuple(a):Fraction(c,r['denominator']) for a,c in r['terms']}
                  for r in records}
        assert {k:d for k,d in total.items() if d}==expected
        print(core,'cyclic coefficients',sorted(expected),'agree with',last-first+1,'individual graphs',flush=True)
    assert all_indices==list(range(3,50))


if __name__=='__main__':
    check()
