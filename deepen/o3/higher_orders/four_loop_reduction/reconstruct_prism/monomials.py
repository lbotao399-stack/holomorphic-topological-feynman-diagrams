"""Reusable exact integrals for every prism monomial in the arbitrary-source kernels."""
import argparse
from collections import defaultdict
from fractions import Fraction as F
from functools import cache
from math import factorial,prod
from pathlib import Path
import json
import re
import subprocess
from tempfile import TemporaryDirectory

import sympy as sp
from evaluate import HERE,BUBBLES,MASTER,Rational,auxiliary,bubble_product,parse_coefficients
from regenerate import edge_permutations

BATCH=200

def form_name(first,cut):
    return f'monomials{first:04}c{cut}'

def make_inputs(record,destination):
    edges=record['edges']
    perms=edge_permutations(tuple(map(tuple,edges)))
    for cut in (0,1):
        order=[cut]+[e for e in range(9) if e!=cut]
        vertices=[[] for _ in range(6)]
        for j,e in enumerate(order):
            u,v=edges[e]
            label='Q' if j==0 else f'p{j}'
            vertices[u].append(label)
            vertices[v].append('-'+label)
        vertex='*'.join('vx('+','.join(v)+')' for v in vertices)
        for first in range(0,len(record['support']),BATCH):
            name=form_name(first,cut)
            terms=[]
            for key,power in enumerate(record['support'][first:first+BATCH],first):
                a=max(tuple(power[p[i]] for i in order) for p in perms)
                prop='*'.join(f'1/(p{i}.p{i})^{n+1}' for i,n in enumerate(a) if i)
                terms.append(f'+{prod(factorial(n) for n in a)}*tag^{key+2}*{prop}')
            text='#-\n#include- forcer.h\nS tag;\nL F='+vertex+'*(\n'+'\n'.join(terms)+'\n);\n#call Forcer()\n'
            (destination/(name+'.frm')).write_text(text+f'#write <{name}.out> "%E",F\n.end\n')


def raw_vectors(text):
    text=''.join(text.split()).replace('\\','')
    pattern=r'([^+]*?\*tag\^(\d+)\*rat\([^()]*?\))(?=\+|$)'
    found=re.findall(pattern,text)
    assert '+'.join(v for v,_ in found)==text
    terms=defaultdict(list)
    for term,tag in found:
        terms[int(tag)-2].append(term.replace(f'*tag^{tag}',''))
    return {key:parse_coefficients('+'.join(vals)) for key,vals in terms.items()}


@cache
def regularized_masters(auxiliary_cuts):
    # Laurent coefficients of the bubble products and of the master value selected
    # by an auxiliary cut difference. Keep each log and zeta term until cancellation.
    delta=auxiliary(*auxiliary_cuts,HERE/'masters')
    table={}
    replacement=defaultdict(int)
    master_upper=1
    for name,js in BUBBLES.items():
        low,factor,taylor=bubble_product(js)
        series={low+i:factor*coef for i,coef in enumerate(taylor)}
        table[name]=series
        multiplier=Rational(0)-delta[name]/delta[MASTER]
        assert multiplier.valuation()+max(series)>=master_upper
        for power,v in multiplier.series(master_upper-low).items():
            for j,t in series.items():
                if power+j<=master_upper:
                    replacement[power+j]+=v*t
    table[MASTER]=dict(replacement)
    symbols=sp.symbols('P E L')
    return {key:{i:{t:F(c) for t,c in sp.Poly(sp.expand(v).xreplace(
                       {sp.pi:symbols[0],sp.EulerGamma:symbols[1],sp.log(2):symbols[2]}),*symbols).terms() if c}
                     for i,v in vals.items()} for key,vals in table.items()}


def verify_auxiliary():
    names=list(BUBBLES)
    expected=[sp.Rational(94208,1715),-sp.Rational(258505869452,4990490505),
              (-138764493638084616019+31022124567494860800*sp.pi**2)/121015527116108625]
    for cuts in ((0,1),(0,2),(1,2)):
        table=regularized_masters(cuts)
        def expression(name):
            return {k:sum(sp.Rational(c)*sp.pi**a*sp.EulerGamma**b*sp.log(2)**d
                          for (a,b,d),c in terms.items()) for k,terms in table[name].items()}
        top,bottom=expression(MASTER),expression(names[1])
        assert bottom[-1]==sp.pi**2/6994593273507840
        answer=[]
        for j in range(3):
            v=top.get(j-1,0)-sum(bottom.get(k,0)*answer[j-1-k] for k in range(j))
            answer.append(sp.factor(v/bottom[-1]))
        assert all(sp.expand(a-b)==0 for a,b in zip(answer,expected))


def evaluate_vector(vector,auxiliary_cuts):
    table=regularized_masters(auxiliary_cuts)
    if MASTER in vector:
        delta=auxiliary(*auxiliary_cuts,HERE/'masters')
        assert vector[MASTER].valuation()>=delta[MASTER].valuation()
    result=defaultdict(F)
    for name,c in vector.items():
        master=table[name]
        assert c.valuation()+max(master)>=0
        for j,f in c.series(-min(master)).items():
            for i,coefs in master.items():
                if i+j>0:continue
                for monomial,v in coefs.items():
                    result[i+j,monomial]+=F(f)*v
    keep=[(0,(2,0,0)),(0,(4,0,0))]
    assert all(k in keep or v==0 for k,v in result.items())
    return tuple(result[k] for k in keep)


def monomial_values(source,master_folder):
    periods=[]
    for first in range(0,len(source['support']),BATCH):
        batches=[]
        for cut,aux in [(0,(0,1)),(1,(1,2))]:
            name=form_name(first,cut)
            vectors=raw_vectors((master_folder/(name+'.out')).read_text())
            ids=range(first,min(first+BATCH,len(source['support'])))
            assert set(vectors)==set(ids)
            batches.append([evaluate_vector(vectors[i],aux) for i in ids])
        assert batches[0]==batches[1]
        periods.extend(batches[0])
        print('monomials',len(periods),flush=True)
    return periods


def contract(source,periods):
    return [[[sum(F(c,source['denominator'])*periods[i][j] for i,c in column) for j in range(2)]
             for column in g] for g in source['graphs']]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('work',type=Path)
    parser.add_argument('--generate-input',action='store_true')
    parser.add_argument('--save',type=Path)
    parser.add_argument('--form',type=Path)
    parser.add_argument('--forcer',type=Path)
    args=parser.parse_args()
    source=json.loads(args.source.read_text())
    if args.generate_input:
        make_inputs(source,args.work)
        return
    if not args.save:
        parser.error('--save is required for reconstructed periods')
    if bool(args.form)!=bool(args.forcer):
        parser.error('--form and --forcer are supplied together')
    if args.form:
        args.work.mkdir(parents=True,exist_ok=True)
        with TemporaryDirectory(prefix='prism-monomials-') as temp:
            folder=Path(temp)
            for name in ('forcer.h','forcer'):
                (folder/name).symlink_to(args.forcer.resolve()/name)
            make_inputs(source,folder)
            for first in range(0,len(source['support']),BATCH):
                for cut in (0,1):
                    name=form_name(first,cut)
                    with (folder/'progress.log').open('w') as log:
                        job=subprocess.run([str(args.form.resolve()),'-q',name+'.frm'],cwd=folder,stdout=log,stderr=subprocess.STDOUT)
                    if job.returncode:
                        raise RuntimeError((folder/'progress.log').read_text()[-2000:])
                    (args.work/(name+'.out')).write_bytes((folder/(name+'.out')).read_bytes())
    verify_auxiliary()
    p=monomial_values(source,args.work)
    kernel=contract(source,p)
    # Machin's identity, with alternating-series rational enclosures for both arctangents.
    def atan_interval(q):
        t=sum(F((-1)**k,(2*k+1)*q**(2*k+1)) for k in range(50))
        return t,t+F(1,101*q**101)
    l5,h5=atan_interval(5);l239,h239=atan_interval(239)
    pi_low,pi_high=16*l5-4*h239,16*h5-4*l239
    assert all(a+b*(pi_low if b>=0 else pi_high)**2>0 for a,b in p)
    bank={'edges':source['edges'],'pi_powers':[2,4],
          'periods':[[u,*map(str,v)] for u,v in zip(source['support'],p)]}
    (HERE/'monomial_periods.json').write_text(json.dumps(bank,separators=(',',':'))+'\n')
    args.save.parent.mkdir(parents=True,exist_ok=True)
    graph_data={'start':source['start'],'powers':source['powers'],
                'graphs':[[list(map(str,coef)) for coef in g] for g in kernel]}
    args.save.write_text(json.dumps(graph_data,separators=(',',':'))+'\n')


if __name__=='__main__':
    main()
