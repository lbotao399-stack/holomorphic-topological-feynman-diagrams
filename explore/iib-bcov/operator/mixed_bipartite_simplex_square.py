"""Integrate the raw current square by collecting its marked core-length monomials."""
from array import array
from collections import defaultdict
from ctypes import CDLL,c_int,c_size_t,c_uint64,c_void_p
from functools import cache
from math import factorial,lcm,prod
from pathlib import Path
from time import perf_counter
import gzip
import os

from flint import fmpq,fmpz_mpoly_ctx
from mixed_current_basis import adapted_roots
from mixed_prism_rung_sum import ColorMoments,PrimitiveChamber,sparse_flow
from mixed_prism_one_star import METRIC,PRIME,residue
from simplex_gram_moments import gram_table,native,pack,simplex_variables,sparse_rows


def by_metric(poly,paths,dimensions):
    occupied=[i for i,path in enumerate(paths) if len(path)>1]
    inside=[i for h in occupied for i in paths[h]]
    ctx=fmpz_mpoly_ctx.get(tuple('x'+str(i) for i in range(len(inside))))
    substitution=simplex_variables(dimensions)
    data=defaultdict(dict)
    for ex,c in poly.terms():
        degrees=[sum(int(ex[k]) for k in path) for path in paths]
        if max(degrees)>=32:raise OverflowError('core length exceeds one exponent lane')
        key=sum(n<<(5*i) for i,n in enumerate(degrees))
        data[key][tuple(int(ex[i]) for i in inside)]=int(c)
    metrics=list(data)
    powers,rows=sparse_rows(ctx.from_dict(data[k]).compose(*substitution) for k in metrics)
    return metrics,powers,rows


@cache
def kernel():
    path=os.environ['MIXED_SIMPLEX_GRAM_LIB']
    library=CDLL(path)
    square,fold,read,free=(library.simplex_square,library.simplex_bipartite_fold,
                          library.simplex_terms,library.free_simplex_terms)
    square.argtypes=[c_void_p,c_void_p,c_int,c_void_p,c_void_p,c_void_p,c_void_p,c_void_p,c_int,c_int]
    fold.argtypes=[c_void_p,c_uint64,c_void_p,c_void_p,c_void_p,c_void_p]
    square.restype=fold.restype=c_void_p
    read.argtypes=[c_void_p,c_void_p,c_void_p];read.restype=c_size_t
    free.argtypes=[c_void_p];free.restype=None
    return square,fold,read,free


def fold_for(order):
    if len(order)==9:return kernel()[1]
    if len(order)!=108:raise ValueError('invalid core orbit action')
    lib=CDLL(os.environ['MIXED_SIMPLEX_GRAM_LIB'])
    fold=lib.simplex_prism_fold
    fold.argtypes=[c_void_p,c_uint64,c_void_p,c_void_p,c_void_p,c_void_p]
    fold.restype=c_void_p
    return fold


def check_prism_powers(data):
    from itertools import combinations
    import networkx as nx
    import numpy as np
    from four_loop import PRISM
    limits=[]
    for size in range(2,9):
        for subset in combinations(range(9),size):
            graph=nx.Graph(PRISM[e] for e in subset)
            h=len(subset)-len(graph)+1
            if h>0 and nx.is_connected(graph) and not nx.has_bridges(graph):
                limits.append((subset,13*h))
    shifts=np.arange(32,-1,-4,dtype=np.uint64)
    for start in range(0,len(data),16384):
        keys=np.fromiter((k for k,v in data[start:start+16384]),dtype=np.uint64)
        powers=(keys[:,None]>>shifts)&15
        for subset,threshold in limits:
            if np.any(powers[:,subset].sum(axis=1)<=threshold):
                raise ArithmeticError('a prism weight fails strict graph-subset convergence')


def root_moments(root,weighted,paths,dimensions,order,metrics):
    hk,hx,hard=by_metric(root,paths,dimensions)
    sk,sx,soft=by_metric(weighted,paths,dimensions)
    if any(2*max((key>>(5*i))&31 for key in hk)+max((key>>(5*i))&31 for key in sk)
            +len(paths[i])-1>=32 for i in range(9)):
        raise OverflowError('the combined core exponent exceeds one lane')
    soft_identity=[((i,1),) for i in range(len(sx))]
    gram=gram_table(hx,hard,sx,soft_identity,dimensions)
    if gram is None:raise OverflowError('the simplex Gram exceeds the exact native budget')
    den,table=gram
    _,_,free_gram=native()
    try:
        bound=den*sum(abs(c) for row in hard for j,c in row)**2*sum(abs(c) for row in soft for j,c in row)
        if bound.bit_length()>126:raise OverflowError('the total numerator exceeds 126 bits')
        square,_,read,free=kernel()
        fold=fold_for(order)
        arrays=(array('Q',hk),array('Q',sk),*pack(soft))
        h,s,off,idx,coeff=(x.buffer_info()[0] for x in arrays)
        offset=sum((len(path)-1)<<(5*i) for i,path in enumerate(paths))
        native_order=array('i',order)
        metric=array('i',[n for row in metrics for n in row])
        totals=defaultdict(int)
        checks=array('i',[0,0])
        pending=[(0,len(hk))]
        pair_count=0
        metric_count=0
        while pending:
            first,last=pending.pop()
            sizes=array('Q',[0,0])
            result=square(table,h,len(sk),s,off,idx,coeff,sizes.buffer_info()[0],first,last)
            if not result:
                middle=(first+last)//2
                if middle==first:raise MemoryError('one hard metric row exceeds the exact block limit')
                pending.extend(((middle,last),(first,middle)))
                continue
            try:
                pair_count+=sizes[0]
                metric_count+=sizes[1]
                orbit_count=array('Q',[0])
                folded=fold(result,offset,native_order.buffer_info()[0],metric.buffer_info()[0],
                            checks.buffer_info()[0],orbit_count.buffer_info()[0])
            finally:
                free(result)
            try:
                keys=array('Q',[0])*orbit_count[0]
                vals=array('Q',[0])*(2*orbit_count[0])
                n=read(folded,keys.buffer_info()[0],vals.buffer_info()[0])
            finally:
                free(folded)
            for i in range(n):
                hi=vals[2*i+1]
                totals[keys[i]]+=vals[2*i]+((hi-(1<<64) if hi>>63 else hi)<<64)
        return den,totals,checks,(len(hk),len(sk),pair_count,metric_count,len(totals))
    finally:
        free_gram(table)


def project(marks,color,destination,core_index=4):
    from mixed_bipartite_matching import matrix_orderings
    os.environ.setdefault('MIXED_LAPLACE_DET','1')
    probe=PrimitiveChamber(core_index,marks,color)
    _,probe.sheared,probe.flow=sparse_flow(probe.graph,color,probe.sheared)
    batch=ColorMoments(marks,(),core_index,True)
    paths=batch.data[5]
    t=batch.data[1].gens()
    voltage=sum(k*prod(t[i] for i in ex) for ex,k in probe.voltage_terms)
    groups=[]
    for root,soft in adapted_roots(probe,batch.data):
        for group in groups:
            if root==group[0]:
                group[1]+=soft
                break
        else:groups.append([root,soft])
    if core_index==4:
        order=matrix_orderings(probe.core)[0]
    elif core_index==3:
        from four_loop import PERMS
        from quartic_planar_periods import reduction
        moves,name,keep=reduction(probe.core)
        assert not moves and name=='prism' and len(PERMS)==12
        order=tuple(keep[k] for transform in PERMS for k in transform)
    else:raise ValueError('the primitive core must be the prism or K3,3')
    metrics=(METRIC,(11,5,2,7,9,6,8,3,4))
    original=PrimitiveChamber(core_index,marks,color)
    expected=[original.numerator(metric,check=True) for metric in metrics]
    actual=[0,0]
    common=1
    full=defaultdict(int)
    def individual():
        for i,(root,soft) in enumerate(groups):
            start=perf_counter()
            result=root_moments(root,soft*voltage,paths,batch.dimensions,order,metrics)
            print('Metric square',i+1,'of',len(groups),result[3],'seconds',round(perf_counter()-start,2),flush=True)
            yield result
    available=hasattr(CDLL(os.environ['MIXED_SIMPLEX_GRAM_LIB']),'fit_metric_shifts')
    if os.environ.get('MIXED_JOIN_VOLTAGE','1' if available else '0')=='1':
        from mixed_metric_voltage import joined_moments
        chunks=joined_moments(groups,voltage,paths,batch.dimensions,order,metrics)
    else:chunks=individual()
    for den,bank,checks,stats in chunks:
        new=lcm(common,den)
        if new!=common:
            for key in full:full[key]*=new//common
            common=new
        for key,value in bank.items():full[key]+=value*(common//den)
        for k,c in enumerate(checks):actual[k]=(actual[k]+c*pow(den,-1,PRIME))%PRIME
    factor=fmpq(str(probe.weight))/2
    assert [x*residue(factor)%PRIME for x in actual]==expected
    scale=factor/(common*factorial(12))
    data=[]
    for key,value in full.items():
        if value:
            powers=tuple(1+((key>>(5*j))&31) for j in reversed(range(9)))
            assert sum(powers)==52 and max(powers)<13
            if core_index==4:
                a=[sum(powers[j:j+3]) for j in (0,3,6)]
                b=[sum(powers[j::3]) for j in range(3)]
                assert max(a+b)<26 and all(52-a[i]-b[j]+powers[3*i+j]>13 for i in range(3) for j in range(3))
            packed=sum(n<<(4*j) for j,n in enumerate(reversed(powers)))
            data.append((packed,value*prod(factorial(n-1) for n in powers)))
    if core_index==3:check_prism_powers(data)
    destination=Path(destination)
    pending=destination.with_suffix('.partial')
    with gzip.open(pending,'wt',compresslevel=1) as out:
        out.writelines(f'{packed} {value*scale}\n' for packed,value in sorted(data))
    pending.replace(destination)
    print('Nonzero calibrated core orbits',len(data),flush=True)


def publish(word,result):
    import fcntl
    import json
    from mixed_bare_symbolic import classified_orbits
    from mixed_bipartite_matching_periods import columns,package

    marks=tuple(map(int,word))
    if len(set(marks))<4:return
    counts={''.join(map(str,c)):n for m,c,n,_ in classified_orbits(4)[0] if m==marks}
    rows=result['colors']
    assert counts=={tag:r['nonzero_hard_orders'] for tag,r in rows.items()}
    period=tuple(sum(p) for p in zip(*(columns(r) for r in rows.values()),strict=True))
    assert result['total']==dict(package(period),nonzero_hard_orders=sum(counts.values()),
                                 weights_sha256=result['total']['weights_sha256'])
    out=Path(__file__).with_name('mixed_bipartite_higher_periods.json')
    fd=os.open(out.parent,os.O_RDONLY)
    try:
        fcntl.flock(fd,fcntl.LOCK_EX)
        all_rows=json.loads(out.read_text()) if out.exists() else {}
        if word in all_rows:
            assert all_rows[word]==result
            return
        all_rows[word]=result
        pending=out.with_suffix('.partial')
        pending.write_text(json.dumps(all_rows,sort_keys=True,indent=2)+'\n')
        pending.replace(out)
    finally:
        os.close(fd)


def full_placement(word,folder):
    import json
    import sys
    from subprocess import run
    from mixed_bare_symbolic import classified_orbits
    from mixed_bipartite_matching import collect_colors
    from mixed_bipartite_matching_periods import columns,digest,package
    from mixed_bipartite_star_bank import measure

    marks=tuple(map(int,word))
    selected=[(color,count) for m,color,count,_ in classified_orbits(4)[0] if m==marks]
    if not selected:raise ValueError('the placement has no surviving canonical color')
    folder=Path(folder)
    folder.mkdir(parents=True,exist_ok=True)
    output=folder/(word+'.json')
    result=json.loads(output.read_text()) if output.exists() else {'colors':{}}
    def save():
        pending=output.with_suffix('.partial')
        pending.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
        pending.replace(output)
    cache_dir=os.environ.get('MIXED_POSITIVE_SCALAR_CACHE')
    if cache_dir is None:
        root=next((p for p in (folder,*folder.parents) if p.name=='.three-path-integrals'),None)
        if root is not None:cache_dir=root/'scalar-functional'
    for c,count in selected:
        tag=''.join(map(str,c))
        src=folder/f'{word}-{tag}.gz'
        if not src.exists():
            run([sys.executable,'-Bu',__file__,word,tag,str(src)],check=True)
        signature=digest(src)
        entry=result['colors'].get(tag)
        if entry is None:
            if cache_dir:
                from mixed_bipartite_scalar_functional import measure as integer_dot
                period=integer_dot(cache_dir,src)
            else:period=measure(folder/'positive',13,folder,src,10000)
            result['colors'][tag]=dict(package(period),weights_sha256=signature,nonzero_hard_orders=count)
            save()
            print('Integrated color',word,tag,*period,flush=True)
        else:
            assert entry==dict(package(columns(entry)),weights_sha256=signature,nonzero_hard_orders=count)
    combined=folder/(word+'-direct.gz')
    collect_colors(word,combined)
    direct=measure(folder/'positive',13,folder,combined,10000)
    added=tuple(sum(p) for p in zip(*(columns(r) for r in result['colors'].values()),strict=True))
    assert tuple(map(str,added))==direct
    result['total']=dict(package(direct),weights_sha256=digest(combined),
                         nonzero_hard_orders=sum(n for c,n in selected))
    save()
    publish(word,result)
    print('Complete placement',word,result['total'],flush=True)


if __name__=='__main__':
    from argparse import ArgumentParser
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--core',choices=('bipartite','prism'),default='bipartite')
    parser.add_argument('word',help='six sorted core path labels of a canonical placement')
    parser.add_argument('color_or_folder',help='one color word or the scratch folder for the full placement')
    parser.add_argument('destination',nargs='?',help='output coefficient bank for a single color')
    args=parser.parse_args()
    from tempfile import TemporaryDirectory
    from subprocess import run
    import shlex
    with TemporaryDirectory() as folder:
        if not os.environ.get('MIXED_SIMPLEX_GRAM_LIB'):
            lib=Path(folder)/'simplex-square.so'
            run([*shlex.split(os.environ.get('CXX','c++')),'-O3','-std=c++17','-shared','-fPIC',
                 str(Path(__file__).with_name('simplex_gram_moments.cpp')),'-o',str(lib)],check=True)
            os.environ['MIXED_SIMPLEX_GRAM_LIB']=str(lib)
        if args.destination is None:
            if args.core=='prism':parser.error('give a color word and coefficient-bank destination for the prism')
            full_placement(args.word,args.color_or_folder)
        else:project(tuple(map(int,args.word)),tuple(map(int,args.color_or_folder)),args.destination,
                     3 if args.core=='prism' else 4)
