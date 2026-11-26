"""Cyclic second disk representative on two divergence-free type-I primary polyvectors."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from functools import cache


@cache
def wedge(arguments):
    """Expand exterior products of sums of input momenta."""
    result = Counter()
    for choices in product(*arguments):
        labels = [label for label, _ in choices]
        if len(set(labels)) != len(labels):
            continue
        sign = (-1)**sum(i > j for pos,i in enumerate(labels) for j in labels[pos+1:])
        for _, factor in choices:
            sign *= factor
        result[tuple(sorted(labels))] += sign
    return {key:value for key,value in result.items() if value}


# External gauge labels are nonnegative. -2 and -1 are the two closed momenta.
CLOSED=(-2,-1)
def monomial(arguments, left, right, scale, out):
    def select(i):
        return ((i,1),) if i in CLOSED else arguments[i]
    for a,x in wedge(tuple(select(i) for i in left)).items():
        for b,y in wedge(tuple(select(i) for i in right)).items():
            if -2 not in a and -1 not in b:  # divergence-free external tensors
                out[a,b] += scale*x*y


def hochschild(column, arity):
    labels = [((i,1),) for i in range(arity+1)]
    left,right = column
    out = Counter()
    monomial(labels[1:],left,right,1,out)
    for i in range(arity):
        args = labels[:i]+[tuple(labels[i]+labels[i+1])]+labels[i+2:]
        monomial(args,left,right,(-1)**(i+1),out)
    monomial(labels[:-1],left,right,(-1)**(arity+1),out)
    return out


def cyclic(column, arity, moving):
    labels = [((i,1),) for i in range(arity)]
    spectator = tuple((i,-1) for i in [*range(arity),*moving])
    out = Counter()
    monomial(labels,*column,1,out)
    monomial(labels[1:]+[spectator],*column,-(-1)**arity,out)
    return out


def reflection(column, arity):
    # Reversing all skew matrices contributes (-1)^(arity+1); permuting the
    # arity suspended open slots after the spectator contributes its Koszul sign.
    labels=[((i,1),) for i in range(arity)]
    sign=(-1)**(arity+1+arity*(arity-1)//2)
    out=Counter()
    monomial(labels,*column,1,out)
    monomial(labels[::-1],*column,-sign,out)
    return out


def bracket(k,l,moving):
    # Plane-wave coefficients for HKR(k) composed with HKR(l), minus the graded reverse.
    total = k+l-1
    labels = [((i,1),) for i in range(total)]
    out = Counter()
    for outer,inner,reverse in [(k,l,False),(l,k,True)]:
        inner_momentum=-2 if reverse else -1
        for i in range(outer):
            inner_slots=tuple(range(i,i+inner))
            fused=tuple((j,1) for j in inner_slots)
            if inner_momentum in moving:
                fused+=((inner_momentum,1),)
            collapsed=labels[:i]+[fused]+labels[i+inner:]
            sign=(-1)**((inner-1)*i)
            if reverse:
                sign *= -(-1)**((k-1)*(l-1))
            for slots,c in wedge(tuple(collapsed)).items():
                pair=(inner_slots,slots) if reverse else (slots,inner_slots)
                out[pair]+=c*sign
    return {key:value for key,value in out.items() if value}


def difference(k,l,moving):
    # The alternating multiderivation symbol is HKR of the Schouten bracket.
    # Only a derivative on one of the two coefficients contributes to this symbol.
    # Ordered combinations enumerate the shuffles; within each block antisymmetry
    # has already supplied the signed permutation of all its covectors.
    from math import comb
    total=k+l-1
    result=Counter(bracket(k,l,moving))
    for (a,b),value in bracket(k,l,moving).items():
        if not (-1 in a or -2 in b):
            continue
        first=[i for i in a if i>=0]
        second=[i for i in b if i>=0]
        original=first+second
        def parity(labels):
            return (-1)**sum(i>j for pos,i in enumerate(labels) for j in labels[pos+1:])
        for subset in combinations(range(total),len(first)):
            complement=tuple(i for i in range(total) if i not in subset)
            left=((-1,) if -1 in a else ())+subset
            right=((-2,) if -2 in b else ())+complement
            sign=parity(original)*parity(subset+complement)
            result[left,right]-=Fraction(sign*value,comb(total,len(subset)))
    return {key:value for key,value in result.items() if value}


def solve(k,l,moving=(),free=()):
    arity=k+l-2
    labels=tuple(range(arity))
    left_labels=((-1,) if -1 in moving else ())+labels
    right_labels=((-2,) if -2 in moving else ())+labels
    columns=[(a,b) for a in combinations(left_labels,k) for b in combinations(right_labels,l)
             if (set(a)|set(b)).issuperset(labels)]
    rows={}
    def equation(key):
        return rows.setdefault(key, [0]*(len(columns)+1))
    target=difference(k,l,moving)
    for i,column in enumerate(columns):
        maps=[('Hochschild',hochschild(column,arity)),('rotation',cyclic(column,arity,moving)),
              ('reflection',reflection(column,arity))]
        for condition, mapping in maps:
            for key,value in mapping.items():
                if value:
                    equation((condition,key))[i]=value
    for key,value in target.items():
        equation(('Hochschild',key))[-1]=value
    # Canonical free coefficients zero, exact echelon elimination on these integer rows.
    basis={}
    for row in rows.values():
        row=[Fraction(x) for x in row]
        for pivot,known in sorted(basis.items()):
            if row[pivot]:
                value=row[pivot]
                row=[a-value*b for a,b in zip(row,known)]
        pivot=next((i for i,x in enumerate(row[:-1]) if x),None)
        if pivot is None:
            assert not row[-1]
        else:
            value=row[pivot]
            basis[pivot]=[x/value for x in row]
    solution=[Fraction(0)]*len(columns)
    free_indices=[i for i in range(len(columns)) if i not in basis]
    if free:
        if len(free)!=len(free_indices):
            raise ValueError('one coefficient per free column is required')
        for i,value in zip(free_indices,free):
            solution[i]=Fraction(value)
    for pivot,row in sorted(basis.items(),reverse=True):
        solution[pivot]=row[-1]-sum(x*y for x,y in zip(row[pivot+1:-1],solution[pivot+1:]))
    # Replay every monomial of the full symbols, including cyclicity.
    for row in rows.values():
        assert sum(x*y for x,y in zip(row,solution))==row[-1]
    return len(columns),len(basis),{c:v for c,v in zip(columns,solution) if v}


if __name__=='__main__':
    for pair,moving in [((1,1),CLOSED),((1,3),CLOSED),((1,5),(-2,)),((3,3),()),
                        ((3,5),()),((3,3),CLOSED),((3,5),(-2,))]:
        size,rank,solution=solve(*pair,moving)
        print('ranks',pair,'moving',moving,'candidate columns',size,'linear rank',rank,'representative')
        print('nonzero terms',len(solution))
    # At dimension five both constant rank-five tensors are proportional.
