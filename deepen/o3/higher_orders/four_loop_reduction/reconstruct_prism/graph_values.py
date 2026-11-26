"""Tagged FORM reductions of all individual prism graphs, using fixed auxiliary cut differences."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from functools import cache
from pathlib import Path
import gzip
import json
import re
import subprocess
import tempfile

import networkx as nx
import sympy as sp
from evaluate import HERE, parse_coefficients, period_from_vector, prism_coefficients
from regenerate import edge_permutations, form_input

DATA=HERE.parent/'reconstruct_cyclic'
GROUPS=[range(i,min(i+4,39)) for i in range(11,39,4)]


def records(group):
    for i in group:
        for row in json.loads((DATA/f'graph{i:02}.json').read_text()):
            yield i,row['N_power'],row


def full_record(record):
    return {**record,'h':[0]*9,'terms':[(tuple(a+b for a,b in zip(t,record.get('h',[0]*9),strict=True)),c) for t,c in record['terms']]}


def batch_input(group,cut):
    parts=[]
    prefix=None
    for idx,k,row in records(group):
        name,raw=form_input(full_record(row),cut)
        header,body=raw.split('*(\n')
        monomials=body.partition('\n);')[0]
        if prefix is None:
            prefix=header
        assert prefix==header
        parts.append(f'+tag^{100*idx+k}*('+monomials+')')
    name=f'graphs{group.start:02}c{cut}'
    prefix=prefix.replace('L F=','S tag;\nL F=')
    text=prefix+'*(\n'+'\n'.join(parts)+'\n);\n#call Forcer()\n'
    return name,text+f'#write <{name}.out> "%E",F\n.end\n'


def vector_strings(text):
    text=''.join(text.split()).replace('\\','')
    # The plus signs in rational coefficients are enclosed in rat(...).
    parts=re.findall(r'([^+]*?\*tag\^(\d+)\*rat\([^()]*?\))(?=\+|$)',text)
    assert '+'.join(term for term,_ in parts)==text
    out=defaultdict(list)
    for term,tag in parts:
        out[int(tag)].append(term.replace(f'*tag^{tag}',''))
    return {divmod(tag,100):'+'.join(terms) for tag,terms in out.items()}


def graph_coefficients(master_folder=HERE/'graphs'):
    values={}
    for group in GROUPS:
        for cut,aux in [(0,(0,1)),(1,(1,2))]:
            f=master_folder/f'graphs{group.start:02}c{cut}.out.gz'
            rows=vector_strings(gzip.decompress(f.read_bytes()).decode())
            assert set(rows)=={(i,k) for i,k,_ in records(group)}
            for key,text in sorted(rows.items()):
                answer=period_from_vector(parse_coefficients(text),aux)
                if cut:
                    assert sp.expand(answer-values[key])==0,key
                else:
                    values[key]=answer
    for k,total in prism_coefficients().items():
        assert sp.expand(sum(v for (_,power),v in values.items() if power==k)-total)==0
    return values


def run_form(directory,engine,name,text):
    (directory/(name+'.frm')).write_text(text)
    with (directory/(name+'.log')).open('w') as log:
        call=subprocess.run([str(engine),'-q',name+'.frm'],cwd=directory,stdout=log,stderr=subprocess.STDOUT)
    if call.returncode:
        raise RuntimeError((directory/(name+'.log')).read_text()[-2000:])
    return (directory/(name+'.out')).read_text()


@cache
def stored_periods():
    data=json.loads((HERE/'monomial_periods.json').read_text())
    return tuple(map(tuple,data['edges'])),{tuple(e):(Fraction(a),Fraction(b)) for e,a,b in data['periods']}


@cache
def proper_subgraphs(edges):
    result=[]
    for mask in range(1,2**len(edges)-1):
        selected=tuple(i for i in range(len(edges)) if mask>>i&1)
        g=nx.Graph()
        g.add_nodes_from(range(6))
        g.add_edges_from(edges[i] for i in selected)
        loops=len(selected)-6+nx.number_connected_components(g)
        if loops:
            result.append((selected,loops))
    return result


def scalar_period(record,form=None,forcer=None):
    """Contract a prism numerator with the exact period bank. FORCER handles unseen orbits."""
    if (form is None)!=(forcer is None):
        raise ValueError('supply FORM and FORCER together')
    prism_edges,bank=stored_periods()
    if tuple(map(tuple,record['edges']))!=prism_edges:
        raise ValueError('expected the prism in source edge order')
    perms=edge_permutations(prism_edges)
    numerator=Counter()
    for exponents,c in full_record(record)['terms']:
        if len(exponents)!=9 or any(not isinstance(i,int) or i<0 for i in exponents) or sum(exponents)!=13:
            raise ValueError('expected nine nonnegative integer exponents of total degree thirteen')
        key=min(tuple(exponents[i] for i in p) for p in perms)
        numerator[key]+=Fraction(c,record['denominator'])
    pi2=pi4=Fraction()
    missing=[]
    for key,c in numerator.items():
        if not c:
            continue
        if key in bank:
            a,b=bank[key]
            pi2+=c*a
            pi4+=c*b
        else:
            if any(2*sum(key[i]+1 for i in subset)<=11*loops for subset,loops in proper_subgraphs(prism_edges)):
                raise ValueError(f'divergent monomial {key}')
            missing.append((key,c))
    known=pi2*sp.pi**2+pi4*sp.pi**4
    if not missing:
        return known
    if form is None:
        raise ValueError(f'{len(missing)} monomial orbits require FORM and FORCER')
    with tempfile.TemporaryDirectory(prefix='prism-period-') as temp:
        d=Path(temp)
        for name in ('forcer.h','forcer'):
            (d/name).symlink_to(Path(forcer).resolve()/name)
        row={**record,'terms':missing,'denominator':1,'h':[0]*9}
        results=[]
        for cut,aux in [(0,(0,1)),(1,(1,2))]:
            name,raw=form_input({**row,'N_power':0},cut)
            text=run_form(d,Path(form).resolve(),name,raw)
            results.append(period_from_vector(parse_coefficients(text),aux))
        assert sp.expand(results[0]-results[1])==0
        return sp.expand(known+results[0])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate-input',type=Path)
    parser.add_argument('--master-folder',type=Path,default=HERE/'graphs')
    args=parser.parse_args()
    if args.generate_input:
        for group in GROUPS:
            for cut in (0,1):
                name,raw=batch_input(group,cut)
                (args.generate_input/(name+'.frm')).write_text(raw)
        return
    values=graph_coefficients(args.master_folder)
    saved=json.loads((HERE/'graph_values.json').read_text())
    for i in range(11,39):
        by_rank={row['N_power']:row for row in json.loads((DATA/f'graph{i:02}.json').read_text())}
        for k,(a,b) in zip(saved['powers_N'],saved['graphs'][i-11],strict=True):
            contracted=scalar_period(by_rank[k]) if k in by_rank else 0
            certified=sp.Rational(a)*sp.pi**2+sp.Rational(b)*sp.pi**4
            assert sp.expand(certified-contracted)==0
            assert sp.expand(certified-values.get((i,k),0))==0
        print(i,' '.join(str(values.get((i,k),0)) for k in (0,2,4)))


if __name__=='__main__':
    main()
