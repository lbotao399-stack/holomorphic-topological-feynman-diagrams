"""Recompute the nine stored prism master vectors with FORM and FORCER."""
from argparse import ArgumentParser
from collections import defaultdict
from fractions import Fraction
from functools import cache
from math import factorial,prod
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import subprocess

import networkx as nx
from evaluate import HERE,prism_coefficients


@cache
def edge_permutations(edges):
    positions={frozenset(e):i for i,e in enumerate(edges)}
    prism=nx.Graph(edges)
    automorphisms=nx.algorithms.isomorphism.GraphMatcher(prism,prism).isomorphisms_iter()
    return tuple(tuple(positions[frozenset((p[u],p[v]))] for u,v in edges) for p in automorphisms)


def form_input(record,strategy):
    edges=record['edges']
    permutations=edge_permutations(tuple(map(tuple,edges)))
    cut=1 if strategy==1 else 0
    edge_order=[cut]+[i for i in range(9) if i!=cut]
    half_edges=[[] for _ in range(6)]
    for i,k in enumerate(edge_order):
        u,v=edges[k]
        name='Q' if i==0 else f'p{i}'
        half_edges[u].append(name)
        half_edges[v].append('-'+name)
    vertex='*'.join('vx('+','.join(row)+')' for row in half_edges)
    weights=defaultdict(Fraction)
    for exponents,numerator in record['terms']:
        choose=min if strategy==2 else max
        powers=choose(tuple(exponents[p[i]] for i in edge_order) for p in permutations)
        weights[powers]+=Fraction(numerator,record['denominator'])*prod(factorial(i) for i in exponents)
    parts=[]
    for powers,c in sorted(weights.items()):
        if c:
            props='*'.join(f'1/(p{i}.p{i})^{a+1}' for i,a in enumerate(powers) if i)
            parts.append(f'+({c})*{props}')
    name=f'c{strategy}-n{record["N_power"]}'
    source='#-\n#include- forcer.h\nL F='+vertex+'*(\n'+'\n'.join(parts)+'\n);\n#call Forcer()\n'
    return name,source+f'#write <{name}.out> "%E",F\n.end\n'


def main():
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--form',type=Path,required=True)
    parser.add_argument('--forcer',type=Path,required=True)
    args=parser.parse_args()
    engine=args.form.resolve()
    records=json.loads((HERE.parent/'reconstruct_cyclic/core11-38.json').read_text())
    reference=prism_coefficients()
    with TemporaryDirectory(prefix='prism-cut-') as temp:
        folder=Path(temp)
        for file in ('forcer.h','forcer'):
            (folder/file).symlink_to(args.forcer.resolve()/file)
        for strategy in range(3):
            for record in records:
                name,text=form_input(record,strategy)
                source=folder/(name+'.frm')
                source.write_text(text)
                with (folder/(name+'.log')).open('w') as log:
                    job=subprocess.run([str(engine),'-q',source.name],cwd=folder,stdout=log,stderr=subprocess.STDOUT)
                if job.returncode:
                    raise RuntimeError((folder/(name+'.log')).read_text()[-2000:])
                def normalized_output(directory):
                    return ''.join((directory/(name+'.out')).read_text().split()).replace('\\','')
                if normalized_output(folder)!=normalized_output(HERE/'masters'):
                    raise AssertionError(f'{name} differs from stored rational master vector')
        assert prism_coefficients(folder)==reference
        print(reference)


if __name__=='__main__':
    main()
