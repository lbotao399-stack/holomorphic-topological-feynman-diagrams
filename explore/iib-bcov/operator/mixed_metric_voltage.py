"""Apply the common voltage after combining the soft Gram factors of different current squares."""
from array import array
from collections import defaultdict
from ctypes import CDLL,c_int,c_int64,c_size_t,c_uint64,c_void_p
from functools import cache
from math import lcm,prod
from time import perf_counter
import os

from mixed_bipartite_simplex_square import by_metric,fold_for,kernel
from simplex_gram_moments import gram_table,native,pack


@cache
def functions():
    lib=CDLL(os.environ['MIXED_SIMPLEX_GRAM_LIB'])
    names=('metric_moment_pairs','empty_metric_moments','fit_metric_shifts','add_metric_shifts',
           'rescale_metric_moments','apply_voltage','free_metric_moments')
    args=((c_void_p,c_void_p,c_int,c_int,c_void_p),(c_int,),
          (c_void_p,c_void_p,c_int,c_void_p,c_size_t),(c_void_p,c_void_p,c_int,c_void_p,c_void_p),
          (c_void_p,c_int64),(c_void_p,c_int,c_void_p,c_void_p,c_void_p,c_void_p,c_size_t,c_size_t,c_void_p),
          (c_void_p,))
    result=[]
    for name,sig,out in zip(names,args,(c_void_p,c_void_p,c_int,c_size_t,None,c_void_p,None),strict=True):
        f=getattr(lib,name);f.argtypes=list(sig);f.restype=out;result.append(f)
    return result


def core_soft(poly,paths):
    extra=[j for path in paths for j in path[1:]]
    lengths=[sum(poly.context().gens()[j] for j in path) for path in paths]
    terms={}
    rebuilt=poly.context().constant(0)
    for ex,c in poly.terms():
        if any(ex[j] for j in extra):continue
        degree=[int(ex[path[0]]) for path in paths]
        assert sum(degree)==2
        terms[sum(n<<(5*i) for i,n in enumerate(degree))]=int(c)
        rebuilt+=c*prod(x**n for x,n in zip(lengths,degree) if n)
    assert rebuilt==poly
    return terms


def joined_moments(groups,voltage,paths,dimensions,order,metrics,voltage_data=None,prepared=None):
    pairs,empty,fit,add,rebase,apply_voltage,free_bundle=functions()
    _,_,read,free=kernel()
    fold=fold_for(order)
    _,_,free_gram=native()
    vk,sx,volts=voltage_data if voltage_data is not None else by_metric(voltage,paths,dimensions)
    q=len(sx)
    cap=18000000//q
    vs=(array('Q',vk),*pack(volts))
    va=[a.buffer_info()[0] for a in vs]
    jac=sum((len(path)-1)<<(5*i) for i,path in enumerate(paths))
    order=array('i',order)
    metric=array('i',[n for row in metrics for n in row])
    norm=sum(abs(c) for row in volts for j,c in row)
    identity=[((i,1),) for i in range(q)]
    acc=empty(q)
    accum_count=0
    common=1
    total_bound=0
    def flush():
        nonlocal acc,accum_count,common,total_bound
        todo=[(0,accum_count)]
        bank=defaultdict(int);checks=array('i',[0,0]);terms=0
        while todo:
            i,j=todo.pop()
            count=array('Q',[0])
            raw=apply_voltage(acc,len(vk),*va,i,j,count.buffer_info()[0])
            if not raw:
                mid=(i+j)//2
                if mid==i:raise MemoryError('voltage row exceeds the core-monomial limit')
                todo.extend(((mid,j),(i,mid)))
                continue
            try:
                terms+=count[0]
                res=fold(raw,jac,order.buffer_info()[0],metric.buffer_info()[0],checks.buffer_info()[0],
                         count.buffer_info()[0])
            finally:free(raw)
            try:
                keys=array('Q',[0])*count[0]
                vals=array('Q',[0])*(2*count[0])
                n=read(res,keys.buffer_info()[0],vals.buffer_info()[0])
            finally:free(res)
            for k in range(n):
                hi=vals[2*k+1]
                bank[keys[k]]+=vals[2*k]+((hi-(1<<64) if hi>>63 else hi)<<64)
        den=common;size=accum_count
        free_bundle(acc);acc=empty(q)
        accum_count=0;common=1;total_bound=0
        return den,bank,checks,(size,terms,len(bank))
    try:
        for r,(root,soft) in enumerate(groups):
            start=perf_counter()
            hk,hx,hard=prepared if prepared is not None else by_metric(root,paths,dimensions)
            shift=core_soft(soft,paths)
            if any(2*max((key>>(5*i))&31 for key in hk)+max((key>>(5*i))&31 for key in vk)+
                   max((key>>(5*i))&31 for key in shift)+len(paths[i])-1>=32 for i in range(9)):
                raise OverflowError('combined length exponent exceeds one lane')
            width=max(1,12000000//len(hx)**2)
            if q>width:
                for start_col in range(0,q,width):
                    stop=min(q,start_col+width)
                    subrows=[[(k-start_col,c) for k,c in v if start_col<=k<stop] for v in volts]
                    data=vk,sx[start_col:stop],subrows
                    yield from joined_moments([(root,soft)],voltage,paths,dimensions,order,metrics,
                                              data,(hk,hx,hard))
                continue
            grams=gram_table(hx,hard,sx,identity,dimensions)
            if grams is None:raise OverflowError(f'exact Gram bound exceeded at {len(hx)} hard monomials and {q} soft columns')
            den,table=grams
            try:
                hs=array('Q',hk);shifts=array('Q',shift);weights=tuple(shift.values())
                bound=den*sum(abs(c) for row in hard for j,c in row)**2*sum(map(abs,weights))*norm
                if bound.bit_length()>126:raise OverflowError('the lifted root exceeds 126 bits')
                todo=[(0,len(hk))]
                while todo:
                    i,j=todo.pop();count=array('Q',[0])
                    piece=pairs(table,hs.buffer_info()[0],i,j,count.buffer_info()[0])
                    if not piece:
                        mid=(i+j)//2
                        if mid==i:raise MemoryError('the paired moments exceed the native limit')
                        todo.extend(((mid,j),(i,mid)))
                        continue
                    try:
                        new=lcm(common,den)
                        if accum_count and (not fit(acc,piece,len(shift),shifts.buffer_info()[0],cap) or
                                (total_bound*(new//common)+bound*(new//den)).bit_length()>126 or
                                (new//common).bit_length()>62 or (new//den*max(map(abs,weights))).bit_length()>62):
                            yield flush()
                            new=den
                        if not fit(acc,piece,len(shift),shifts.buffer_info()[0],cap):
                            mid=(i+j)//2
                            if mid==i:raise MemoryError('one shifted pair block exceeds the native limit')
                            todo.extend(((mid,j),(i,mid)))
                            continue
                        if accum_count and new!=common:rebase(acc,new//common)
                        w=array('q',[c*(new//den) for c in weights])
                        accum_count=add(acc,piece,len(shift),shifts.buffer_info()[0],w.buffer_info()[0])
                        total_bound=total_bound*(new//common)+bound*(new//den)
                        common=new
                    finally:free_bundle(piece)
            finally:free_gram(table)
            print('Joined soft metric',r+1,'of',len(groups),len(hk),accum_count,
                  'seconds',round(perf_counter()-start,2),flush=True)
        if accum_count:yield flush()
    finally:free_bundle(acc)
