"""Symbolic harmonic-number representation of a complete bipartite angular shell."""
from collections import defaultdict
from functools import cache
from itertools import product
from math import comb,factorial

from flint import fmpq as Q
from quartic_bipartite_angular import columns_and_scale
from quartic_bipartite_sums import excluded_sum, rational_triangle


# A term (a,p,w) represents (s+a/2)^(-p) H_s(w), where H is a strict
# finite multiple harmonic sum. Polynomial terms have a=0 and p<=0.
def compact(data):
    return tuple(sorted((a,p,w,c) for (a,p,w),c in data.items() if c))


def add(*terms):
    result=defaultdict(Q)
    for factor,value in terms:
        for a,p,w,c in value:result[a,p,w]+=factor*c
    return compact(result)


def rational(a,p):
    if p<0:
        return tuple((0,-j,(),Q(comb(-p,j))*Q(a,2)**(-p-j)) for j in range(-p+1) if a or j==-p)
    return ((a if p else 0,p,(),Q(1)),)


def constant(c):
    return ((0,0,(),Q(c)),) if c else ()


@cache
def product_poles(a,p,b,q):
    if p>q:return product_poles(b,q,a,p)
    if q<=0:return rational(0,p+q)
    if p<=0:return add(*[(comb(-p,j)*Q(-b,2)**(-p-j),rational(b,q-j)) for j in range(1-p)])
    if a==b:return rational(a,p+q)
    result=[]
    for u,v,k,l in ((a,b,p,q),(b,a,q,p)):
        d=Q(v-u,2)
        for j in range(k):result.append(((-1)**j*comb(l+j-1,j)*d**(-l-j),rational(u,k-j)))
    return add(*result)


@cache
def stuffle(u,v):
    if not u:return (v,)
    if not v:return (u,)
    return (tuple((u[0],)+t for t in stuffle(u[1:],v))
            +tuple((v[0],)+t for t in stuffle(u,v[1:]))
            +tuple((u[0]+v[0],)+t for t in stuffle(u[1:],v[1:])))


def multiply(left,right):
    result=defaultdict(Q)
    for a,p,u,x in left:
        for b,q,v,y in right:
            for k,l,_,z in product_poles(a,p,b,q):
                weight=x*y*z
                for w in stuffle(u,v):result[k,l,w]+=weight
    return compact(result)


@cache
def shift_word(w,offset):
    result=((0,0,w,Q(1)),)
    if not w or not offset:return result
    if offset>0:
        indices=range(1,offset+1);sign=1
    else:
        indices=range(offset+1,1);sign=-1
    return add((1,result),*[(sign,multiply(rational(2*j,w[0]),shift_word(w[1:],j-1))) for j in indices])


@cache
def power_primitive(k):
    # Forward difference of the resulting polynomial is s^k.
    coefficients=[Q(0)]*(k+2)
    for j in range(k,-1,-1):
        lower=sum((coefficients[i]*comb(i,j) for i in range(j+2,k+2)),Q(0))
        coefficients[j+1]=(int(j==k)-lower)/(j+1)
    return compact({(0,-j,()):c for j,c in enumerate(coefficients)})


def shift_rational(poly,k):
    return add(*[(c,rational(a+2*k,p)) for a,p,w,c in poly])


@cache
def primitive_term(a,p,w):
    if a%2:raise ValueError('half-integer pole in integer primitive')
    if p>0:
        head=shift_word((p,)+w,a//2-1)
        if not w:return head
        difference=add((1,shift_word(w,0)),(-1,shift_word(w,a//2-1)))
        remainder=multiply(rational(a,p),difference)
    else:
        poly=power_primitive(-p)
        head=multiply(poly,shift_word(w,0))
        if not w:return head
        remainder=add((-1,multiply(multiply(shift_rational(poly,1),rational(2,w[0])),shift_word(w[1:],0))))
    return add((1,head),*[(v,primitive_term(b,q,word)) for b,q,word,v in remainder])


@cache
def hsum(n,w):
    if not w:return Q(1)
    if n<=0:return Q(0)
    return hsum(n-1,w)+hsum(n-1,w[1:])/Q(n)**w[0]


def evaluate(expression,s):
    return sum(c*Q(2*s+a,2)**(-p)*hsum(s,w) for a,p,w,c in expression)


@cache
def interval(p,a):
    if p>0:return add((1,shift_word((p,),a)),(-1,constant(hsum(a-1,(p,)))))
    poly=power_primitive(-p)
    return add((1,shift_rational(poly,a+1)),(-1,constant(evaluate(poly,a))))


@cache
def diagonal(p,q,a,b):
    if p<=0:
        return add(*[(comb(-p,j)*(-1)**j,multiply(rational(2*(a+b),p+j),interval(q-j,b)))
                     for j in range(1-p)])
    if q<=0:return diagonal(q,p,b,a)
    return add(*[(comb(p+q-h-1,other-1),multiply(rational(2*(a+b),p+q-h),interval(h,offset)))
                 for own,other,offset in ((p,q,a),(q,p,b)) for h in range(1,own+1)])


@cache
def pair(p,q,a,b):
    if (q,b)<(p,a):return pair(q,p,b,a)
    full=multiply(interval(p,a),interval(q,b))
    anti=add(*[(c,primitive_term(k,l,w)) for k,l,w,c in diagonal(p,q,a,b)])
    base=max(a,b)+2
    offset=Q(str(excluded_sum(p,q,a,b,base)))-evaluate(anti,base)
    return add((1,full),(-1,anti),(-1,constant(offset)))


@cache
def triple(axes):
    if axes!=tuple(sorted(axes)):return triple(tuple(sorted(axes)))
    (p,a),(q,b),(r,c)=axes
    d=a+b+c
    if p==0:return pair(q,r,b,c)
    if p<0:
        return add(*[(comb(-p,j)*comb(-p-j,k)*(-1)**(j+k)*2**(-p-j-k),
                      multiply(rational(d,p+j+k),pair(q-j,r-k,b,c)))
                     for j in range(1-p) for k in range(1-p-j)])
    return multiply(add((Q(1,2),triple(((p-1,a),(q,b),(r,c)))),
                        (Q(1,2),triple(((p,a),(q-1,b),(r,c)))),
                        (Q(1,2),triple(((p,a),(q,b),(r-1,c))))),rational(d,1))


@cache
def directional_binomial(q,a):
    # Expand binom(s-n+q-2,q-2) in n+a, retaining polynomials in s.
    out=[constant(Q(1,factorial(q-2)))]
    for j in range(1,q-1):
        out=[add((1,multiply(v,rational(2*(a+j),-1))),(-1,previous))
             for v,previous in zip(out+[()], [()]+out)]
    return tuple(out)


def shell_expression(p,q,partition=0):
    columns,scale=columns_and_scale(p,q,partition)
    lam=q-1
    final=defaultdict(Q)
    for exterior in (False,True):
        dim=sum(b[2] for b in columns)-q if exterior else q
        factors=[]
        for b in columns:
            entries={}
            for l,h,a,e,c in rational_triangle(2*q-sum(b) if exterior else b[0],b[1],b[2],q):
                for j,weight in enumerate(directional_binomial(q,a)):
                    key=(l,h,a,e-j)
                    entries[key]=add((1,entries.get(key,())),(Q(str(c)),weight))
            factors.append([(*k,v) for k,v in entries.items() if v])
        groups=defaultdict(lambda:())
        for first,second,third in product(*factors):
            rows=(first,second,third)
            l,h=dim+sum(i[0] for i in rows),sum(i[1] for i in rows)
            axes=tuple(sorted((i[3],i[2]) for i in rows))
            product_weight=multiply(multiply(first[4],second[4]),third[4])
            value=multiply(product_weight,rational(2*l,h+1))
            key=axes
            groups[key]=add((1,groups[key]),(Q((-1)**h*factorial(h),2**(h+1)),value))
        for axes,weight in groups.items():
            for a,k,w,c in multiply(weight,triple(axes)):final[a,k,w]+=c
    out=compact(final)
    for d in range(lam+1,2*lam):out=multiply(out,rational(2*d,-1))
    norm=Q(str(scale))*Q(2*lam*factorial(lam-1),factorial(lam)*factorial(2*lam-1))
    return add((norm,out),)


def controls():
    import json
    from pathlib import Path
    from quartic_bipartite_sums import shell,triple_moment
    for axes in (((2,1),(1,3),(3,2)),((0,2),(1,1),(-2,2)),((-2,3),(-1,4),(1,2))):
        value=triple(axes)
        for n in (3,7,12):assert evaluate(value,n)==Q(str(triple_moment(n,axes)))
    source=json.loads(Path(__file__).with_name('quartic_terminal_periods.json').read_text())
    for row in source['K33']:
        p=row['powers']
        value=shell_expression(p,3)
        for n in (2,5):assert evaluate(value,n)==Q(str(shell(p,3,n)))
        print('symbolic shell terms',len(value),'half-integer',sum(a%2 for a,k,w,c in value),
              'words',sorted(set(w for a,k,w,c in value)),flush=True)


if __name__=='__main__':
    controls()
