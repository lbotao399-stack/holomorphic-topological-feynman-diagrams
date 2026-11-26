"""Exact common-denominator Kronecker contraction of homogeneous exterior path polynomials."""
from collections import Counter
from functools import cache
from itertools import groupby
from math import lcm,prod

from flint import fmpq,fmpz,fmpz_mpoly_ctx


class ExteriorProduct:
    def __init__(self,polynomials,occupied,context):
        self.context=context
        self.variables=[i for i in range(context.nvars()) if i not in occupied]
        self.terms=[]
        self.high=[]
        self.low=[]
        self.norm=[]
        self.degree=[]
        for poly in polynomials:
            rows=[]
            for ex,c in poly.terms():
                assert c.denominator==1
                powers=tuple(int(ex[i]) for i in self.variables)
                rows.append((powers,int(c)))
            assert rows and len({sum(ex) for ex,_ in rows})==1
            self.terms.append(rows)
            self.high.append(tuple(map(max,zip(*(ex for ex,_ in rows)))))
            self.low.append(tuple(map(min,zip(*(ex for ex,_ in rows)))))
            self.degree.append(sum(rows[0][0]))
            self.norm.append(sum(abs(c) for _,c in rows))

    def contract(self,terms,mode):
        active=[(f,c) for (f,pole),c in terms.items() if c]
        if not active:return self.context.constant(0)
        degree=sum(self.degree[i] for i in active[0][0])
        denominator=lcm(*(int(c.denominator) for _,c in active))
        weighted=[(f,int(c*denominator)) for f,c in active]
        high=[0]*len(self.variables)
        low=[degree]*len(self.variables)
        bound=0
        for (i,j,k),c in weighted:
            assert self.degree[i]+self.degree[j]+self.degree[k]==degree
            bound+=abs(c)*self.norm[i]*self.norm[j]*self.norm[k]
            high=[max(m,a+b+d) for m,a,b,d in zip(high,self.high[i],self.high[j],self.high[k])]
            low=[min(m,a+b+d) for m,a,b,d in zip(low,self.low[i],self.low[j],self.low[k])]
        sizes=[b-a+1 for a,b in zip(low,high)]
        drop=max(range(len(sizes)),key=sizes.__getitem__)
        keep=[i for i in range(len(sizes)) if i!=drop]
        strides=[];cells=1
        for i in keep:strides.append(cells);cells*=sizes[i]
        frequency=Counter(i for f,c in weighted for i in set(f))
        ordered=sorted((tuple(sorted(f,key=lambda i:(-frequency[i],i))),c) for f,c in weighted)
        def product_sum(factors,zero):
            value=zero
            for i,g in groupby(ordered,key=lambda item:item[0][0]):
                mid=zero
                for j,h in groupby(g,key=lambda item:item[0][1]):
                    last=sum((c*factors(f[2]) for f,c in h),zero)
                    mid+=factors(j)*last
                value+=factors(i)*mid
            return value
        width=(bound.bit_length()+1+7)//8
        if mode=='integer' or cells*width>200000:
            compact=fmpz_mpoly_ctx.get(tuple('v'+str(k) for k in keep) or ('constant',))
            @cache
            def reduced(i):
                return compact.from_dict({tuple(ex[k] for k in keep) or (0,):c for ex,c in self.terms[i]})
            value=product_sum(reduced,compact.constant(0))
            result={}
            for powers,c in value.terms():
                ex=[0]*self.context.nvars()
                for k,v in zip(keep,powers):ex[self.variables[k]]=int(v)
                ex[self.variables[drop]]=degree-sum(ex)
                assert min(ex)>=0
                result[tuple(ex)]=fmpq(int(c),denominator)
            return self.context.from_dict(result)
        bits=8*width
        @cache
        def packed(i):
            return sum((fmpz(c)<<(bits*sum(ex[k]*stride for k,stride in zip(keep,strides)))
                        for ex,c in self.terms[i]),fmpz(0))
        value=product_sum(packed,fmpz(0))
        value>>=bits*sum(low[k]*stride for k,stride in zip(keep,strides))
        if not value:return self.context.constant(0)
        # The L1 bound makes each coefficient a unique balanced radix digit.
        sign=1 if value>0 else -1
        absolute=int(abs(value))
        data=absolute.to_bytes((absolute.bit_length()+7)//8,'little')
        base=1<<bits
        result={}
        carry=0
        for n,start in enumerate(range(0,len(data),width)):
            digit=int.from_bytes(data[start:start+width],'little')+carry
            carry=int(2*digit>=base)
            digit-=carry*base
            if not digit:continue
            assert abs(digit)<=bound
            ex=[0]*self.context.nvars()
            for k,stride in zip(keep,strides):
                ex[self.variables[k]]=low[k]+n//stride%sizes[k]
            ex[self.variables[drop]]=degree-sum(ex)
            assert min(ex)>=0
            result[tuple(ex)]=fmpq(sign*digit,denominator)
        assert not carry
        return self.context.from_dict(result)
