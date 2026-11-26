"""Certified ordinary multiple-zeta evaluation of an integrated angular shell."""
from collections import Counter,defaultdict
from functools import cache
from itertools import combinations
from math import comb,factorial

from flint import fmpq as Q
from quartic_harmonic_shell import add,evaluate,primitive_term,shell_expression,stuffle


@cache
def regularized(word):
    # Constant term at H_s(1)=0 in the finite harmonic stuffle algebra.
    if not word or word[0]>1:return {word:Q(1)}
    terms=Counter(stuffle((1,),word[1:]))
    lead=terms.pop(word)
    out=defaultdict(Q)
    for other,c in terms.items():
        for key,v in regularized(other).items():out[key]-=Q(c,lead)*v
    return dict(out)


@cache
def shift_series(j,degree,logs):
    """Coefficient of x^degree in (1-x)^(-j) log(1-x)^logs."""
    if logs:
        return -sum((shift_series(j,degree-i,logs-1)/i for i in range(1,degree+1)),Q(0))
    return Q(comb(j+degree-1,degree) if j else int(degree==0))


def shifted(coefficient,j,degree):
    return [(h-k,w,c*comb(h,k)*shift_series(j,degree,k))
            for (h,w),c in coefficient.items() for k in range(min(h,degree)+1)]


@cache
def asymptotic(word,last):
    # Expand H_s(word) in s^{-1}, with L=log(s)+EulerGamma. The
    # coefficient at level j obeys (d/dL-j)P_j from its backward difference.
    if not word:return [{(0,()):Q(1)}]+[{} for _ in range(last)]
    output=[]
    low=asymptotic(word[1:],max(0,last+1-word[0]))
    for j in range(last+1):
        rhs=defaultdict(Q)
        previous=j+1-word[0]
        if previous>=0:
            for i in range(previous+1):
                for h,w,c in shifted(low[i],i,previous-i):rhs[h,w]+=c
        for i,old in enumerate(output):
            for h,w,c in shifted(old,i,j+1-i):rhs[h,w]+=c
        level=defaultdict(Q)
        if j==0:
            level.update({(0,w):c for w,c in regularized(word).items()})
            for (h,w),c in rhs.items():level[h+1,w]+=c/(h+1)
        else:
            for (h,w),c in rhs.items():
                for t in range(h+1):
                    level[h-t,w]-=Q(factorial(h),factorial(h-t)*j**(1+t))*c
        output.append({key:v for key,v in level.items() if v})
    return output


@cache
def words(weight):
    if weight==0:return ((),)
    return tuple((i,)+tail for i in range(1,weight+1) for tail in words(weight-i))


def binary(word):
    return tuple(k for n in word for k in (0,)*(n-1)+(1,))


def decode(bits):
    result=[];n=0
    for k in bits:
        n+=1
        if k:result.append(n);n=0
    assert not n
    return tuple(result)


@cache
def shuffles(first,second):
    if not first:return (second,)
    if not second:return (first,)
    return (tuple((first[0],)+v for v in shuffles(first[1:],second))
            +tuple((second[0],)+v for v in shuffles(first,second[1:])))


def relations(weight):
    # Iterated-integral reversal, sum formula and shuffle/stuffle products.
    domain=[w for w in words(weight) if w[0]>1]
    for w in domain:
        other=decode(tuple(1-k for k in reversed(binary(w))))
        yield [(w,1),(other,-1)]
    for depth in range(2,weight):
        yield [(w,1) for w in domain if len(w)==depth]+[((weight,),-1)]
    for k in range(1,weight//2+1):
        for left in words(k):
            if left[0]==1 and left!=(1,):continue
            for right in words(weight-k):
                if right[0]==1:continue
                identity=Counter(stuffle(left,right))
                identity.subtract(decode(v) for v in shuffles(binary(left),binary(right)))
                assert all(not c or w[0]>1 for w,c in identity.items())
                yield identity.items()


@cache
def relation_basis(weight):
    # Prefer the simple zeta as the uneliminated coefficient in each weight.
    indices=sorted((w for w in words(weight) if w[0]>1),key=lambda w:(-len(w),w))
    order={w:i for i,w in enumerate(indices)}
    basis={}
    for raw in relations(weight):
        coeff=defaultdict(Q)
        for w,c in raw:coeff[w]+=c
        col={w:c for w,c in coeff.items() if c}
        while col:
            key=min(col,key=order.get)
            if key not in basis:
                basis[key]={w:c/col[key] for w,c in col.items()};break
            scale=col[key]
            for w,c in basis[key].items():
                v=col.get(w,Q(0))-scale*c
                if v:col[w]=v
                elif w in col:del col[w]
    return tuple((w,basis[w]) for w in indices if w in basis)


def reduce_words(input_values):
    values=defaultdict(Q,input_values)
    for weight in sorted({sum(w) for w,c in values.items() if w and c}):
        for pivot,row in relation_basis(weight):
            scale=values[pivot]
            if scale:
                for w,c in row.items():values[w]-=scale*c
    return {w:c for w,c in values.items() if c}


def period(p,q,partition=0):
    from quartic_bipartite_angular import columns_and_scale,shell_term
    expression=shell_expression(p,q,partition)
    if any(a%2 for a,h,w,c in expression):raise ValueError('A half-integer shell pole survives.')
    anti=add(*[(c,primitive_term(a,h,w)) for a,h,w,c in expression])
    lower=max([0]+[1-a//2 for a,h,w,c in anti if h>0 and a<=0])
    values=defaultdict(Q)
    for a,h,w,c in anti:
        if h>0:continue
        for (logs,word),v in asymptotic(w,-h)[-h].items():
            if not logs:values[word]+=c*v
    columns,scale=columns_and_scale(p,q,partition)
    for s in range(lower):
        values[()]+=Q(str(scale*sum(shell_term(i,j,s-i-j,columns,q)
                                     for i in range(s+1) for j in range(s-i+1))))
    values[()]-=evaluate(anti,lower)
    return reduce_words(values)


def controls():
    import json
    from pathlib import Path
    rows=json.loads(Path(__file__).with_name('quartic_terminal_periods.json').read_text())['K33']
    for row in rows:
        expected={w:2*Q(c) for w,c in zip(((),(3,),(5,)),row['1,zeta3,zeta5']) if Q(c)}
        for partition in (0,1):
            actual=period(row['powers'],3,partition)
            assert actual==expected,(row,actual,expected)
        print('independent harmonic period',row['powers'],actual,flush=True)


def compare(form,forcer,q=4):
    from quartic_recut_bank import finite_controls
    from math import prod
    for i,row in enumerate(finite_controls(form,forcer,q),1):
        p=row['powers']
        normalization=Q(prod(factorial(a-1) for a in p),factorial(q-1))
        expected={w:Q(c)/normalization for w,c in zip(((),(3,),(5,)),row['1,zeta3,zeta5']) if Q(c)}
        actual=period(p,q)
        assert actual==expected,(p,actual,expected)
        print('orbit',i,p,actual,flush=True)


if __name__=='__main__':
    from argparse import ArgumentParser
    from pathlib import Path
    cli=ArgumentParser(description=__doc__)
    cli.add_argument('--compare',type=Path,nargs=2)
    args=cli.parse_args()
    if args.compare:compare(*args.compare)
    else:controls()
