"""Contract the three exterior factors against the dimension-26 K3,3 period bank.

The projected coefficient carries the reciprocal of the marked-degree factorials.
For each raw monomial the angular kernel carries all nine exponent factorials
divided by 12!, so the marked factors cancel. The three output columns multiply
1, zeta(3), zeta(5), in the same normalization as mixed_bipartite_matching.
"""
from argparse import ArgumentParser
from array import array
from collections import Counter, defaultdict
from contextlib import nullcontext
from ctypes import CDLL, c_int, c_size_t, c_uint64, c_void_p, POINTER
from fractions import Fraction
from functools import cache
from hashlib import sha256
from itertools import combinations, groupby
from math import comb, factorial, lcm, prod
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory
from time import perf_counter
import gc
import gzip
import json
import os
import shlex

import numpy as np
from flint import fmpq, fmpz_mpoly_ctx
from mixed_bipartite_one_star import PrimitiveChamber
from mixed_bipartite_star_bank import ARCHIVE, COMPLETION, verify_archive, verify_member
from mixed_prism_one_star import chamber_rows
from mixed_prism_rung_sum import ColorMoments

HERE=Path(__file__).resolve().parent
FACT=tuple(map(factorial,range(13)))
# Positive powers use the bank's four-bit layout; exterior monomials use
# five-bit lanes so a term outside the convergent range cannot alias another.
# The bipartition has row vertices (0,4,5) and column vertices (1,2,3).
MATRIX=(0,1,2,3,5,7,4,6,8)
NATIVE=r"""
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
#define main original_orbit_driver
#include "quartic_bipartite_orbit_ibp.cpp"
#undef main
#pragma GCC diagnostic pop
extern "C" void chart(size_t n, const uint64_t *raw, uint64_t *out){
    for(size_t i=0;i<n;i++){
        Matrix p{};
        uint64_t r=raw[i];
        for(int j=8;j>=0;j--){p[j]=1+(r&31);r>>=5;}
        out[i]= std::accumulate(p.begin(),p.end(),0)==52 && finite(p,13) ? canonical(p) : 0;
    }
}

#include <queue>
struct DuffyBank {
    static constexpr size_t slots=1<<21;
    std::vector<uint64_t> raw_keys=std::vector<uint64_t>(slots,0);
    std::vector<uint32_t> raw_ids=std::vector<uint32_t>(slots,0);
    std::unordered_map<uint64_t,uint32_t> ids;
    std::vector<uint64_t> keys;
    std::vector<__int128_t> result;
    uint32_t index(uint64_t r){
        uint64_t mix=r;
        mix=(mix^(mix>>30))*0xbf58476d1ce4e5b9ULL;
        mix=(mix^(mix>>27))*0x94d049bb133111ebULL;
        size_t slot=(mix^(mix>>31))&(slots-1);
        if(raw_keys[slot]==r)return raw_ids[slot];
        Matrix a{};
        uint64_t k=r;
        for(int j=8;j>=0;j--){a[j]=1+(k&31);k>>=5;}
        k=std::accumulate(a.begin(),a.end(),0)==52 && finite(a,13) ? canonical(a) : (1ULL<<60)|r;
        auto entry=ids.emplace(k,keys.size());
        if(entry.second){keys.push_back(k);result.push_back(0);}
        raw_keys[slot]=r;
        raw_ids[slot]=entry.first->second;
        return entry.first->second;
    }
};
extern "C" void* duffy_start(){return new DuffyBank;}
extern "C" void duffy_stop(void* bank){delete static_cast<DuffyBank*>(bank);}
struct MetricPair {uint64_t metric; uint32_t left,right;};
extern "C" size_t duffy_square(void* ptr,int nh,int ns,int ha,int sw,int sums,uint64_t offset,
        const uint64_t* hm,const uint64_t* sm,const int* he,const int* hi,const int64_t* hc,
        const int* se,const int* si,const int64_t* sc,const int* power_sum,const int64_t* moments){
    auto& bank=*static_cast<DuffyBank*>(ptr);
    std::fill(bank.result.begin(),bank.result.end(),0);
    // Merge the sorted hard-metric rows with linear storage on wider roots.
    auto later=[](const MetricPair& a,const MetricPair& b){return a.metric>b.metric;};
    std::priority_queue<MetricPair,std::vector<MetricPair>,decltype(later)> heap(later);
    std::vector<MetricPair> pairs;
    size_t pair_count=size_t(nh)*(nh+1)/2;
    bool compact=pair_count>4000000;
    if(compact){
        for(int i=0;i<nh;i++)heap.push({hm[i]+hm[i],uint32_t(i),uint32_t(i)});
    }else{
        pairs.reserve(pair_count);
        for(int i=0;i<nh;i++)for(int j=i;j<nh;j++)pairs.push_back({hm[i]+hm[j],uint32_t(i),uint32_t(j)});
        std::sort(pairs.begin(),pairs.end(),[](const MetricPair& a,const MetricPair& b){return a.metric<b.metric;});
    }
    size_t consumed=0;
    auto front=[&](){return compact?heap.top():pairs[consumed];};
    auto take=[&](){
        auto pair=front();
        if(compact){
            heap.pop();
            auto j=pair.right+1;
            if(j<uint32_t(nh))heap.push({hm[pair.left]+hm[j],pair.left,j});
        }
        ++consumed;
        return pair;
    };
    std::vector<int64_t> square(sums,0);
    std::vector<uint32_t> serial(sums,0),active;
    std::vector<__int128_t> integral(sw);
    uint32_t stamp=0;
    while(consumed<pair_count){
        uint64_t metric=front().metric;
        ++stamp;active.clear();
        do{
            auto pair=take();
            auto i=pair.left,j=pair.right;
            for(int x=he[i];x<he[i+1];x++){
                auto ids=power_sum+int64_t(hi[x])*ha;
                auto c=hc[x]*(1+(i!=j));
                for(int y=he[j];y<he[j+1];y++){
                    int k=ids[hi[y]];
                    if(serial[k]!=stamp){serial[k]=stamp;active.push_back(k);square[k]=0;}
                    square[k]+=c*hc[y];
                }
            }
        }while(consumed<pair_count && front().metric==metric);
        std::fill(integral.begin(),integral.end(),0);
        for(auto k:active)if(square[k]){
            const auto* row=moments+int64_t(k)*sw;
            for(int w=0;w<sw;w++)integral[w]+=__int128_t(square[k])*row[w];
        }
        for(int s=0;s<ns;s++){
            __int128_t value=0;
            for(int t=se[s];t<se[s+1];t++)value+=integral[si[t]]*sc[t];
            if(value){
                auto index=bank.index(offset+metric+sm[s]);
                bank.result[index]+=value;
            }
        }
    }
    return bank.keys.size();
}
extern "C" size_t duffy_read(void* ptr,uint64_t* keys,uint64_t* values){
    auto& b=*static_cast<DuffyBank*>(ptr);
    size_t count=0;
    for(size_t i=0;i<b.keys.size();i++)if(b.result[i]){
        auto v=__uint128_t(b.result[i]);
        keys[count]=b.keys[i];values[2*count]=uint64_t(v);values[2*count+1]=uint64_t(v>>64);++count;
    }
    return count;
}
"""


def native(folder):
    src=folder/'chart.cpp'
    src.write_text(NATIVE)
    lib=folder/'chart.so'
    run([*shlex.split(os.environ.get('CXX','c++')),'-O3','-std=c++17','-shared','-fPIC',
         '-I',str(HERE.parent),str(src),'-o',str(lib)],check=True)
    handle=CDLL(str(lib))
    handle.chart.argtypes=[c_size_t,POINTER(c_uint64),POINTER(c_uint64)]
    handle.chart.restype=None
    handle.duffy_start.argtypes=[]
    handle.duffy_start.restype=c_void_p
    handle.duffy_stop.argtypes=[c_void_p]
    handle.duffy_stop.restype=None
    handle.duffy_square.argtypes=[c_void_p,*([c_int]*5),c_uint64,*([c_void_p]*10)]
    handle.duffy_square.restype=c_size_t
    handle.duffy_read.argtypes=[c_void_p,c_void_p,c_void_p]
    handle.duffy_read.restype=c_size_t
    return handle


def packed(exponents):
    key=0
    for i in MATRIX:
        key=32*key+int(exponents[i])
    return key


class Periods:
    def __init__(self,handle):
        self.handle=handle
        self.canonical={}
        self.values={}

    def charts(self,keys):
        keys=np.ascontiguousarray(keys,dtype=np.uint64)
        values=np.empty_like(keys)
        pointer=POINTER(c_uint64)
        self.handle.chart(keys.size,keys.ctypes.data_as(pointer),values.ctypes.data_as(pointer))
        return values

    def register(self,raw):
        keys=np.fromiter((k for k in raw if k not in self.canonical),dtype=np.uint64)
        if not len(keys):return
        values=self.charts(keys)
        self.canonical.update(zip(keys.tolist(),values.tolist(),strict=True))

    def read(self,keys=None):
        wanted=set(self.canonical.values() if keys is None else keys).difference(self.values,{0})
        if not wanted:return
        verify_archive()
        meta=json.loads((HERE/'mixed_bipartite_star_bank.json').read_text())['completion']
        verify_member(COMPLETION,meta['extra_positive']+meta['extra_contracted'],meta['sha256'])
        for path in (ARCHIVE,COMPLETION):
            with gzip.open(path,'rb') as stream:
                for row in stream:
                    k,tail=row.split(None,1)
                    key=int(k)
                    if key not in wanted:continue
                    wanted.remove(key)
                    self.values[key]=tuple(fmpq(v.decode()) for v in tail.split())
                    if not wanted:break
            if not wanted:break
        if wanted:raise ArithmeticError(f'{len(wanted)} convergent periods missing from the bank')


def integer_kernel(periods):
    denominators=[lcm(*(int(row[i].denominator) for row in periods.values())) for i in range(3)]
    kernel={}
    for k,cols in periods.items():
        fac=prod(FACT[((k>>(4*i))&15)-1] for i in range(9))
        kernel[k]=tuple(int(c*d)*fac for c,d in zip(cols,denominators,strict=True))
    return denominators,kernel


def raw_terms(poly,integral=False):
    terms=[]
    maxima=[0]*9
    for powers,coefficient in poly.terms():
        ex=tuple(map(int,powers))
        maxima=[max(a,b) for a,b in zip(maxima,ex,strict=True)]
        if integral:
            if coefficient.denominator!=1:raise ArithmeticError('Exterior pair has fractional content')
            coefficient=int(coefficient)
        terms.append((packed(ex),coefficient))
    return terms,maxima


class FactorContraction:
    def __init__(self,marked,periods):
        self.marked=marked
        self.periods=periods
        self.products=[]
        self.expanded_pairs=0

    def collect(self,part):
        polynomials=part.polynomials
        @cache
        def pair(i,j):
            self.expanded_pairs+=1
            return raw_terms(polynomials[i]*polynomials[j],integral=True)
        for degree,terms in part.coefficients.items():
            words=[(f,c) for (f,pole),c in terms.items() if c]
            if not words:continue
            frequency=Counter(p for f,_ in words for p in set(combinations(f,2)))
            leaves=defaultdict(list)
            for f,c in words:
                p=max(combinations(f,2),key=lambda x:(frequency[x],-len(polynomials[x[0]])*len(polynomials[x[1]])))
                last=list(f)
                for i in p:last.remove(i)
                leaves[p].append((last[0],c))
            marked=[0]*9
            for path,n in zip(self.marked,degree,strict=True):marked[path]=n
            offset=packed(marked)
            for (i,j),rows in leaves.items():
                linear=sum(c*polynomials[k] for k,c in rows)
                if not linear:continue
                terms3,high3=raw_terms(linear)
                terms2,high2=pair(i,j)
                if any(x+y+d>=32 for x,y,d in zip(high2,high3,marked,strict=True)):
                    raise ArithmeticError('Exterior product exceeds the integer chart layout')
                self.products.append((offset,terms2,terms3,prod(map(factorial,degree))))
                keys={offset+m+n for m,_ in terms2 for n,_ in terms3}
                self.periods.register(keys)

    def value(self):
        self.periods.read()
        total=[fmpq(0) for _ in range(3)]
        invalid=defaultdict(fmpq)
        # Use the common denominator of each angular column for the integer
        # transpose. The marked factorials cancel the simplex volume exactly.
        denom,angular=integer_kernel(self.periods.values)
        kernel={k:angular[v] for k,v in self.periods.canonical.items() if v}
        for offset,pair,leaf,marked_factorial in self.products:
            values=[0,0,0]
            d=lcm(*(int(c.denominator) for n,c in leaf))
            for n,c in leaf:
                shift=offset+n
                left=int(c*d)
                col=[0,0,0]
                for m,b in pair:
                    v=kernel.get(shift+m)
                    if v is None:
                        invalid[shift+m]+=fmpq(left,d*marked_factorial)*b
                        continue
                    for i in range(3):col[i]+=b*v[i]
                for i in range(3):values[i]+=left*col[i]
            for i in range(3):
                total[i]+=fmpq(values[i],d*marked_factorial*FACT[12]*denom[i])
        if any(invalid.values()):raise ArithmeticError('Nonconvergent monomials survive the physical sum')
        return tuple(map(str,total))


class PackedFunctional:
    """Partial Kronecker substitution encodes exterior directions in each integer.

    Homogeneity restores one more direction. Only the other exterior powers
    occur in the symbolic product. The final contraction extracts the middle
    coefficient of a product with each reversed angular functional.
    """
    def __init__(self,marked,periods,packed_coordinates=2):
        self.marked=marked
        self.packed_coordinates=packed_coordinates
        self.free=[i for i in range(9) if i not in marked]
        self.periods=periods
        self.denominator=1
        self.bound=0
        self.maxima=[0]*len(self.free)
        self.slices={}

    def inspect(self,part):
        maxima=[]
        norms=[]
        total_degrees=[]
        for poly in part.polynomials:
            rows=list(poly.terms())
            if any(any(m[i] for i in self.marked) or c.denominator!=1 for m,c in rows):
                raise ArithmeticError('Exterior factors must be unmarked and integral')
            totals={sum(map(int,m)) for m,c in rows}
            if len(totals)!=1:raise ArithmeticError('Inhomogeneous exterior factor')
            total_degrees.append(totals.pop())
            maxima.append([max(int(m[i]) for m,c in rows) for i in self.free])
            norms.append(sum(abs(int(c)) for m,c in rows))
        high=np.array(maxima,dtype=np.int64)
        total_degrees=np.array(total_degrees,dtype=np.int64)
        # A sum of positive L1 bounds covers cancellation between roots as well.
        for degree,terms in part.coefficients.items():
            words=[(f,c) for (f,pole),c in terms.items() if c]
            if not words:continue
            ix=np.array([f for f,c in words],dtype=np.int32)
            if not np.all(total_degrees[ix].sum(axis=1)+sum(degree)==43):
                raise ArithmeticError('Exterior degrees disagree with the projective weight')
            if max(degree)>=32:raise ArithmeticError('Marked degree exceeds the integer chart layout')
            self.maxima=np.maximum(self.maxima,high[ix].sum(axis=1).max(axis=0)).tolist()
            height=max(abs(int(c.numerator)) for f,c in words)
            self.bound+=len(words)*height*max(norms)**3
            self.denominator=lcm(self.denominator,*(int(c.denominator) for f,c in words))

    def accept(self,part):
        old_denominator=self.denominator
        self.inspect(part)
        if max(self.maxima)>=32:
            raise ArithmeticError('Exterior degree exceeds the integer chart layout')
        if not self.slices:
            order=sorted(range(len(self.free)),key=lambda i:self.maxima[i])
            self.drop=self.free[order[-1]]
            self.embed=[self.free[i] for i in order[:self.packed_coordinates]]
            self.variables=[self.free[i] for i in order[self.packed_coordinates:-1]]
            self.ctx=fmpz_mpoly_ctx.get(tuple(f'e{i}' for i in self.variables) or ('constant',))
            self.stride=1
            self.width=0
        stride=self.maxima[self.free.index(self.embed[0])]+1
        width=((self.bound*self.denominator).bit_length()+2+7)//8
        width=self.width if self.width>=width else width+1
        scale=self.denominator//old_denominator
        if self.slices and (width!=self.width or stride!=self.stride):
            # Only the compressed polynomial is recoded. No graph orbit weights
            # are formed. A radix change is needed only when a new root exceeds
            # the accumulated certified support or coefficient bound.
            bits=8*self.width
            radix=1<<bits
            mask=radix-1
            for degree,poly in self.slices.items():
                rows={}
                for ex,c in poly.terms():
                    previous=int(c)
                    encoded=0
                    i=0
                    while previous:
                        digit=previous&mask
                        if digit>=radix//2:digit-=radix
                        previous=(previous-digit)>>bits
                        if digit:
                            j=i%self.stride+stride*(i//self.stride)
                            encoded+=(scale*digit)<<(8*width*j)
                        i+=1
                    rows[tuple(map(int,ex))]=encoded
                self.slices[degree]=self.ctx.from_dict(rows)
        elif self.slices and scale!=1:
            for degree in self.slices:self.slices[degree]*=scale
        self.width=width
        self.stride=stride
        self.length=stride*(1+self.maxima[self.free.index(self.embed[1])]) if len(self.embed)==2 else stride
        self.collect(part)

    def collect(self,part):
        polynomials=[]
        for poly in part.polynomials:
            rows=defaultdict(int)
            for m,c in poly.terms():
                ix=int(m[self.embed[0]])
                if len(self.embed)==2:ix+=self.stride*int(m[self.embed[1]])
                ex=tuple(int(m[i]) for i in self.variables) or (0,)
                rows[ex]+=int(c)<<(8*self.width*ix)
            polynomials.append(self.ctx.from_dict(rows))
        zero=self.ctx.constant(0)
        for degree,terms in part.coefficients.items():
            active=[(f,int(c*self.denominator)) for (f,pole),c in terms.items() if c]
            if not active:continue
            frequency=Counter(i for f,c in active for i in set(f))
            ordered=sorted((tuple(sorted(f,key=lambda i:(-frequency[i],i))),c) for f,c in active)
            poly=self.slices.get(degree,zero)
            for i,first in groupby(ordered,key=lambda item:item[0][0]):
                inner=zero
                for j,group in groupby(first,key=lambda item:item[0][1]):
                    last=sum((c*polynomials[f[2]] for f,c in group),zero)
                    inner+=polynomials[j]*last
                poly+=polynomials[i]*inner
            self.slices[degree]=poly

    def plan(self):
        """Push the angular functional back along the partial substitution."""
        self.records={}
        wanted=set()
        shifts={i:1<<(5*(8-MATRIX.index(i))) for i in range(9)}
        index=np.arange(self.length,dtype=np.int64)
        a=index%self.stride
        b=index//self.stride
        translated=a*shifts[self.embed[0]]
        if len(self.embed)==2:translated+=b*shifts[self.embed[1]]
        self.arguments=0
        for degree,poly in self.slices.items():
            base=[]
            remaining=[]
            encoded=[]
            fixed=[0]*9
            for i,v in zip(self.marked,degree,strict=True):fixed[i]=v
            for ex,value in poly.terms():
                for j,n in zip(self.variables,ex):fixed[j]=int(n)
                base.append(packed(fixed))
                remaining.append(43-sum(fixed))
                encoded.append(int(value))
            if not base:continue
            h=np.array(remaining,dtype=np.int64)[:,None]-a-b
            valid=(h>=0)&(h<=self.maxima[self.free.index(self.drop)])
            raw=np.where(valid,np.array(base,dtype=np.int64)[:,None]+translated+h*shifts[self.drop],0)
            charts=self.periods.charts(raw)
            wanted.update(np.unique(charts).tolist())
            self.arguments+=int(valid.sum())
            self.records[degree]=encoded,charts
        self.slices.clear()
        self.periods.read(wanted)

    def value(self):
        self.plan()
        denom,angular=integer_kernel(self.periods.values)
        # Four signed radix digits give one constant-term pairing. Enlarging
        # the radix makes every intervening product coefficient a balanced digit.
        bits=8*self.width
        wide=bits+max(abs(c).bit_length() for v in angular.values() for c in v)+4
        mask=(1<<wide)-1
        half=1<<(wide-1)
        unit=1<<bits
        def widen_two(n):
            hi=(n+unit//2)>>bits
            return (hi<<wide)+n-(hi<<bits)
        total=[fmpq(0) for _ in range(3)]
        for degree,(encoded,charts) in self.records.items():
            values=[0,0,0]
            for polynomial,keys in zip(encoded,charts.tolist(),strict=True):
                for start in range(0,len(keys),4):
                    chunk=polynomial% (1<<(4*bits))
                    if chunk>=1<<(4*bits-1):chunk-=1<<(4*bits)
                    polynomial=(polynomial-chunk)>>(4*bits)
                    if not chunk:continue
                    hi=(chunk+(1<<(2*bits-1)))>>(2*bits)
                    enlarged=(widen_two(hi)<<(2*wide))+widen_two(chunk-(hi<<(2*bits)))
                    reverse=[0,0,0]
                    for j in range(4):
                        reverse=[r<<wide for r in reverse]
                        k=keys[start+j] if start+j<len(keys) else 0
                        if k:
                            v=angular[k]
                            for col in range(3):reverse[col]+=v[col]
                        else:
                            # On nonconvergent or absent chart entries the
                            # coefficient itself must cancel before integration.
                            r=(enlarged+(1<<(j*wide-1)))>>(j*wide) if j else enlarged
                            if r&mask:raise ArithmeticError('Nonconvergent monomial survives')
                    for col in range(3):
                        center=((enlarged*reverse[col]+(1<<(3*wide-1)))>>(3*wide))&mask
                        values[col]+=center if center<half else center-(1<<wide)
                if polynomial:raise ArithmeticError('Radix support exceeded')
            norm=prod(map(factorial,degree))*FACT[12]*self.denominator
            for i in range(3):total[i]+=fmpq(values[i],norm*denom[i])
        return tuple(map(str,total))


@cache
def duffy_path(ex):
    """Integrate a simplex monomial on independent unit intervals by binomial expansion."""
    value=Fraction(1)
    tail=0
    for degree in reversed(ex):
        value*=sum((Fraction((-1)**j*comb(tail,j),degree+j+1) for j in range(tail+1)),Fraction())
        tail+=1+degree
    return value


def factor_source(part,rows,substitution):
    """Transpose the rank factorization of every marked degree into a sparse metric map."""
    source=defaultdict(lambda:defaultdict(int))
    common=lcm(*(int(row[3].denominator) for row in rows))
    for degree,local,i,c,pole in rows:
        if pole:raise ArithmeticError('Duffy contraction requires an uncut core')
        fixed=[0]*9
        for j,k in zip(part.marked,degree,strict=True):fixed[j]=k
        base=packed(fixed)
        angle=[]
        for m,k in local.compose(*substitution).terms():
            key=0
            for d in m:
                if d>=32:raise OverflowError('angular exponent exceeds a lane')
                key=(key<<5)|int(d)
            angle.append((key,int(k)))
        scale=int(c*common)
        for ex,coefficient in part.polynomials[i].terms():
            if coefficient.denominator!=1:raise ArithmeticError('noninteger primitive exterior polynomial')
            value=scale*int(coefficient)
            row=source[base+packed(ex)]
            for a,h in angle:row[a]+=value*h
    ids={}
    metric=[]
    ends=array('i',[0]);indices=array('i');weights=array('q')
    for k,row in sorted(source.items()):
        used=False
        for a,v in row.items():
            if not v:continue
            coefficient,remainder=divmod(v,common)
            if remainder:raise ArithmeticError('RREF did not reconstruct the integral source')
            indices.append(ids.setdefault(a,len(ids)))
            weights.append(coefficient)
            used=True
        if used:
            metric.append(k)
            ends.append(len(weights))
    return array('Q',metric),tuple(ids),ends,indices,weights


def duffy_moments(hard,soft,dimensions):
    hard=np.array(hard,dtype=np.int64)
    soft=np.array(soft,dtype=np.int64)
    lanes=np.arange(sum(dimensions))*5
    max_h=(hard[:,None]>>lanes&31).max(axis=0)
    max_s=(soft[:,None]>>lanes&31).max(axis=0)
    if np.any(2*max_h+max_s>=32):raise OverflowError('angular exponent sum exceeds a lane')
    adjacent,indices=np.unique(hard[:,None]+hard,return_inverse=True)
    indices=np.ascontiguousarray(indices.reshape(len(hard),len(hard)),dtype=np.int32)
    powers,inverse=np.unique(adjacent[:,None]+soft,return_inverse=True)
    exact=[]
    for key in powers.tolist():
        ex=tuple((key>>i)&31 for i in reversed(range(0,5*sum(dimensions),5)))
        position=0
        value=Fraction(1)
        for size in dimensions:
            value*=duffy_path(ex[position:position+size])
            position+=size
        exact.append(value)
    common=lcm(*(v.denominator for v in exact))
    lifted=np.array([v.numerator*(common//v.denominator) for v in exact],dtype=np.int64)
    moments=np.ascontiguousarray(lifted[inverse].reshape(len(adjacent),len(soft)),dtype=np.int64)
    return common,indices,moments


def angular_at_metric(source,metric):
    from mixed_prism_one_star import PRIME
    times=np.array([metric[i] for i in MATRIX],dtype=np.int64)
    metric_keys,angle,ends,indices,weights=source
    powers=np.array([[(key>>j)&31 for j in reversed(range(0,45,5))] for key in metric_keys],dtype=np.int64)
    factors=np.ones(len(powers),dtype=np.int64)
    for j,t in enumerate(times.tolist()):
        table=np.array([pow(t,k,PRIME) for k in range(32)],dtype=np.int64)
        factors=factors*table[powers[:,j]]%PRIME
    result=np.zeros(len(angle),dtype=np.int64)
    np.add.at(result,np.asarray(indices),np.repeat(factors,np.diff(ends))*(np.asarray(weights)%PRIME))
    return result%PRIME


class DuffyFunctional:
    def __init__(self,part,handle):
        from mixed_prism_one_star import METRIC
        self.part=part
        self.handle=handle
        self.bank=handle.duffy_start()
        self.weights=defaultdict(int)
        self.common=1
        self.metrics=(METRIC,(11,5,2,7,9,6,8,3,4))
        self.checks=[0,0]
        fixed=[0]*9
        for j,d in zip(part.marked,part.dimensions,strict=True):fixed[j]=d
        self.offset=packed(fixed)
        ctx=fmpz_mpoly_ctx.get(tuple(f'u{i}' for i in range(sum(part.dimensions))))
        variables=ctx.gens()
        self.sub=[]
        at=0
        for size in part.dimensions:
            group=variables[at:at+size]
            self.sub.extend((*group,1-sum(group)))
            at+=size

    def root(self,hard,weighted):
        from mixed_prism_one_star import PRIME
        part=self.part
        h=factor_source(part,part.split(hard),self.sub)
        w=factor_source(part,part.split(weighted),self.sub)
        part.ids.clear()
        part.polynomials.clear()
        for i in range(0,45,5):
            peak=2*max((k>>i)&31 for k in h[0])+max((k>>i)&31 for k in w[0])+((self.offset>>i)&31)
            if peak>=32:raise OverflowError('metric degree sum exceeds a lane')
        den,indices,moments=duffy_moments(h[1],w[1],part.dimensions)
        hard_bound=sum(map(abs,h[4]))
        if (2*hard_bound**2).bit_length()>63:raise OverflowError('squared angular coefficients exceed int64')
        if (den*hard_bound**2*sum(map(abs,w[4]))).bit_length()>126:raise OverflowError('period numerator exceeds int128')
        for j,metric in enumerate(self.metrics):
            left=angular_at_metric(h,metric);right=angular_at_metric(w,metric)
            matrix=((moments%PRIME@right)%PRIME)[indices]
            value=int(((left@matrix%PRIME)@left)%PRIME)*pow(den,-1,PRIME)
            jac=prod(metric[k]**d for k,d in zip(part.marked,part.dimensions,strict=True))
            self.checks[j]=(self.checks[j]+jac*value)%PRIME
        arrays=(h[0],w[0],*h[2:],*w[2:])
        count=self.handle.duffy_square(self.bank,len(h[0]),len(w[0]),len(h[1]),len(w[1]),len(moments),
                                      self.offset,*(a.buffer_info()[0] for a in arrays),indices.ctypes.data,moments.ctypes.data)
        keys=array('Q',[0])*count
        values=array('Q',[0])*(2*count)
        rows=self.handle.duffy_read(self.bank,keys.buffer_info()[0],values.buffer_info()[0])
        common=lcm(self.common,den)
        if common!=self.common:
            for k in self.weights:self.weights[k]*=common//self.common
            self.common=common
        factor=common//den
        for i in range(rows):
            hi=values[2*i+1]
            signed=(hi-(1<<64) if hi>>63 else hi)<<64
            self.weights[keys[i]]+=(values[2*i]+signed)*factor
        return len(h[0]),len(w[0]),rows

    def calculate(self,color,progress,reference_weights):
        from mixed_current_basis import adapted_roots
        from mixed_prism_rung_sum import sparse_flow
        from mixed_prism_one_star import PRIME,residue
        part=self.part
        probe=PrimitiveChamber(4,part.marks,color)
        _,probe.sheared,probe.flow=sparse_flow(probe.graph,color,probe.sheared)
        times=part.data[1].gens()
        voltage=sum(k*prod(times[i] for i in ex) for ex,k in probe.voltage_terms)
        groups=[]
        for r,w in adapted_roots(probe,part.data):
            for previous in groups:
                if previous[0]==r:
                    previous[1]+=w
                    break
            else:groups.append([r,w])
        for i,(r,w) in enumerate(groups):
            start=perf_counter()
            shape=self.root(r,w*voltage)
            if progress:print('Duffy root',i+1,'of',len(groups),'metric rows',shape,'seconds',round(perf_counter()-start,2),flush=True)
        scale=fmpq(str(probe.weight))/2
        reference=PrimitiveChamber(4,part.marks,color)
        for value,metric in zip(self.checks,self.metrics,strict=True):
            if value*residue(scale)%PRIME!=reference.numerator(metric,check=True):
                raise ArithmeticError('Independent grounded metric disagrees')
        weights={k:v for k,v in self.weights.items() if v}
        if any(k>>60 for k in weights):raise ArithmeticError('Nonconvergent raw powers survived the full sum')
        bank=Periods(self.handle)
        bank.read(weights)
        den,kernel=integer_kernel(bank.values)
        total=[0,0,0]
        physical_scale=scale/(self.common*FACT[12])
        digest=sha256()
        source=gzip.open(reference_weights,'rb') if reference_weights else nullcontext()
        with source as reference:
            for k in sorted(weights):
                c=weights[k]
                for i in range(3):total[i]+=c*kernel[k][i]
                fac=prod(FACT[((k>>(4*i))&15)-1] for i in range(9))
                row=f'{k} {physical_scale*(c*fac)}\n'.encode()
                if reference is not None and reference.readline()!=row:
                    raise ArithmeticError(f'Independent orbit weight differs at {k}')
                digest.update(row)
            if reference is not None and reference.readline():raise ArithmeticError('Extra reference orbit weights')
        print('Calibrated orbit weights and SHA-256:',len(weights),digest.hexdigest(),flush=True)
        return tuple(str(physical_scale*fmpq(v,d)) for v,d in zip(total,den,strict=True))


def calculate(word,color,handle,method=None,packed_coordinates=1,progress=False,reference_weights=None):
    marks=tuple(map(int,word))
    colors=tuple(map(int,color))
    if len(marks)!=6 or marks!=tuple(sorted(marks)) or any(i not in range(9) for i in marks):
        raise ValueError('A placement requires six sorted core path labels')
    if colors not in {row[0] for row in chamber_rows(marks,4)}:
        raise ValueError('The color is not a canonical chamber of this placement')
    batch=ColorMoments(marks,(),4,adapted=True)
    start=perf_counter()
    if method=='duffy' or (method is None and (reference_weights or len(set(marks))>=3)):
        problem=DuffyFunctional(batch,handle)
        try:
            value=problem.calculate(colors,progress,reference_weights)
        except OverflowError as error:
            if method=='duffy' or reference_weights:raise
            print('Using arbitrary-precision partial substitution:',error,flush=True)
            batch=ColorMoments(marks,(),4,adapted=True)
            method='packed'
        else:
            print('Duffy period columns:',*value,'elapsed seconds:',round(perf_counter()-start,2),flush=True)
            return value
        finally:handle.duffy_stop(problem.bank)
    periods=Periods(handle)
    weight=PrimitiveChamber(4,marks,colors).weight
    use_packed=method=='packed'
    if use_packed:
        problem=PackedFunctional(batch.marked,periods,packed_coordinates)
        roots=0
        def accept(part):
            nonlocal roots
            problem.accept(part)
            roots+=1
            if progress:print('Exterior roots contracted:',roots,'elapsed seconds:',round(perf_counter()-start,2),flush=True)
        batch.add(colors,weight,[],accept)
    else:
        problem=FactorContraction(batch.marked,periods)
        batch.add(colors,weight,[],problem.collect)
    if not batch.nonzero:raise ValueError('The hard color has no nonzero chamber')
    built=perf_counter()
    value=problem.value()
    if not use_packed:print('Distinct exterior factor pairs:',problem.expanded_pairs,flush=True)
    print('Raw arguments and canonical kernels:',problem.arguments if use_packed else len(periods.canonical),len(periods.values),flush=True)
    print('Build and angular contraction seconds:',built-start,perf_counter()-built,flush=True)
    print('Period columns:',*value,flush=True)
    return value


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('word')
    parser.add_argument('color')
    parser.add_argument('--method',choices=('pairs','packed','duffy'))
    parser.add_argument('--packed-coordinates',type=int,choices=(1,2),default=1)
    parser.add_argument('--check',type=Path,help='independent exact period or matching-color JSON')
    parser.add_argument('--progress',action='store_true')
    parser.add_argument('--check-weights',type=Path,help='compare every normalized Duffy weight against its raw archive')
    args=parser.parse_args()
    if args.check_weights and args.method not in (None,'duffy'):parser.error('--check-weights requires the Duffy contraction')
    gc.disable()
    with TemporaryDirectory() as tmp:
        value=calculate(args.word,args.color,native(Path(tmp)),args.method,args.packed_coordinates,args.progress,args.check_weights)
        if args.check:
            reference=json.loads(args.check.read_text())
            if 'period' not in reference:reference=reference[args.word]['colors'][args.color]
            assert list(value)==reference['period']
