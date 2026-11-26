"""Exact Dirichlet inner products for low-degree current polynomials along a core path."""
from collections import defaultdict
from fractions import Fraction as F
from functools import cache
from itertools import combinations_with_replacement
from math import comb, factorial, lcm, prod
from flint import fmpq, fmpq_mpoly_ctx


def plus(x,y):
    return tuple(a+b for a,b in zip(x,y))


@cache
def moment(ex):
    return F(factorial(len(ex)-1)*prod(factorial(i) for i in ex),factorial(len(ex)-1+sum(ex)))


class Simplex:
    def __init__(self,m,degree):
        self.m = m
        self.basis = []
        self.norm = []
        self.indices = []
        # Monomials of the first m-1 fractions, in increasing total degree.
        for total in range(degree+1):
            for indices in combinations_with_replacement(range(m-1),total):
                ex = tuple(indices.count(i) for i in range(m))
                weights = self.expand(ex)
                orthogonal = {ex:F(1)}
                for i,c in weights:
                    for key,v in self.basis[i].items():
                        orthogonal[key] = orthogonal.get(key,0)-c*v
                orthogonal = {k:v for k,v in orthogonal.items() if v}
                self.norm.append(sum(v*moment(plus(k,ex)) for k,v in orthogonal.items()))
                self.indices.append(ex)
                self.basis.append(orthogonal)

    def expand(self,ex):
        return [(i,c) for i,b in enumerate(self.basis) if sum(self.indices[i])<=sum(ex)
                if (c:=sum(v*moment(plus(k,ex)) for k,v in b.items())/self.norm[i])]

    @cache
    def position(self,mode,offset):
        """Multiply an orthogonal basis element by the sum of the first offset fractions."""
        result = defaultdict(F)
        for ex,coef in self.basis[mode].items():
            for i in range(offset):
                new = list(ex)
                new[i] += 1
                for j,c in self.expand(tuple(new)):
                    result[j] += c*coef
        return tuple((j,v) for j,v in result.items() if v)


@cache
def simplex(m,degree):
    return Simplex(m,degree)


@cache
def legendre_product(a,b):
    coeff = [0]*(a+b+1)
    for i in range(a+1):
        for j in range(b+1):
            coeff[i+j] += (-1)**(a+b-i-j)*comb(a,i)*comb(a+i,i)*comb(b,j)*comb(b+j,j)
    return coeff


@cache
def minimum_integral(a,b,c,d):
    return sum(x*y*(F(1,i+2)+F(1,j+2))/(i+j+3)
               for i,x in enumerate(legendre_product(a,b))
               for j,y in enumerate(legendre_product(c,d)))


class Cube:
    """Independent labeled positions, with the same segment polynomial on every order."""
    def __init__(self,m,degree):
        self.indices = [tuple(choice.count(i) for i in range(m-1)) for size in range(degree+1)
                        for choice in combinations_with_replacement(range(m-1),size)]
        self.number = {ex:i for i,ex in enumerate(self.indices)}
        self.norm = [F(1,prod(2*k+1 for k in ex)) for ex in self.indices]

    @cache
    def expand(self,ex):
        # Segment lengths are s_1,s_2-s_1,...,1-s_last at unit path resistance.
        size = len(ex)-1
        polynomial = {(0,)*size:1}
        for j,count in enumerate(ex):
            for _ in range(count):
                new = defaultdict(int)
                for powers,v in polynomial.items():
                    q = list(powers)
                    if j<size:
                        q[j] += 1
                    new[tuple(q)] += v
                    if j:
                        q = list(powers)
                        q[j-1] += 1
                        new[tuple(q)] -= v
                polynomial = new
        out = []
        for i,mode in enumerate(self.indices):
            coeff = sum(c*prod(F(factorial(k)**2,factorial(k-r)*factorial(k+r+1))
                               for k,r in zip(powers,mode))
                        for powers,c in polynomial.items() if all(r<=k for k,r in zip(powers,mode)))
            if coeff:
                out.append((i,coeff/self.norm[i]))
        return tuple(out)

    @cache
    def position(self,mode,offset):
        source = self.indices[mode]
        r = source[offset-1]
        rows = [(mode,F(1,2))]
        for step,coefficient in ((1,r+1),(-1,r)):
            if coefficient:
                target = list(source)
                target[offset-1] += step
                rows.append((self.number[tuple(target)],F(coefficient,2*(2*r+1))))
        return tuple(rows)

    @cache
    def minimum(self,mode,first,second):
        a,b = first-1,second-1
        source = self.indices[mode]
        result = []
        for i,target in enumerate(self.indices):
            if all(target[k]==v for k,v in enumerate(source) if k not in (a,b)):
                v = minimum_integral(source[a],target[a],source[b],target[b])
                if v:
                    result.append((i,v*(2*target[a]+1)*(2*target[b]+1)))
        return tuple(result)


@cache
def cube(m,degree):
    return Cube(m,degree)


def rational(v):
    return fmpq(v.numerator,v.denominator)


class PathMoments:
    def __init__(self,polynomial,paths,degree,all_orders=False,*,ctx=None):
        self.ctx = fmpq_mpoly_ctx.get(tuple(f'T{j}' for j in range(len(paths)))) if ctx is None else ctx
        self.all_orders = all_orders
        self.simplexes = [(cube if all_orders else simplex)(len(p),degree+1) for p in paths]
        entries = {}
        for ex,c in polynomial.terms():
            groups = tuple(tuple(int(ex[j]) for j in path) for path in paths)
            assert all(sum(g)<=degree for g in groups)
            entries[(tuple(map(sum,groups)),groups)] = int(c)
        divisor = 1
        # Apply one finite change-of-basis matrix on each path. Keep integer
        # coefficients until all ten path-position integrations have separated.
        for i,space in enumerate(self.simplexes):
            columns = {modes[i] for _,modes in entries}
            chart = {key:space.expand(key) for key in columns}
            scale = lcm(*(c.denominator for col in chart.values() for _,c in col))
            table = {key:[(j,int(v*scale)) for j,v in col] for key,col in chart.items()}
            new = defaultdict(int)
            for (powers,modes),coef in entries.items():
                for j,v in table[modes[i]]:
                    new[powers,modes[:i]+(j,)+modes[i+1:]] += coef*v
            entries = {key:c for key,c in new.items() if c}
            divisor *= scale
        by_mode = defaultdict(dict)
        for (powers,modes),coef in entries.items():
            by_mode[modes][powers] = fmpq(coef,divisor)
        self.base = {m:self.ctx.from_dict(p) for m,p in by_mode.items()}
        divisor = 1 if all_orders else prod(factorial(len(p)-1) for p in paths)
        self.volume = self.ctx.from_dict({tuple(len(p)-1 for p in paths):fmpq(1,divisor)})
        self.positions = {}
        self.norms = {}

    def position(self,path,offset):
        return self._apply(path,(offset,))

    def minimum(self,path,first,second):
        if not self.all_orders or first==second:
            return self.position(path,min(first,second))
        return self._apply(path,tuple(sorted((first,second))))

    def _apply(self,path,offsets):
        if (path,offsets) not in self.positions:
            new = {}
            space = self.simplexes[path]
            mapping = space.position if len(offsets)==1 else space.minimum
            for modes,p in self.base.items():
                for j,v in mapping(modes[path],*offsets):
                    key = modes[:path]+(j,)+modes[path+1:]
                    new[key] = new.get(key,0)+rational(v)*p
            self.positions[path,offsets] = {key:v for key,v in new.items() if v}
        return self.positions[path,offsets]

    def inner(self,left,right):
        result = self.ctx.constant(0)
        for key in left.keys() & right.keys():
            if key not in self.norms:
                self.norms[key] = rational(prod(s.norm[i] for s,i in zip(self.simplexes,key)))
            result += self.norms[key]*left[key]*right[key]
        return result


def add_scaled(target,source,scale):
    if scale:
        if isinstance(scale,F):
            scale = rational(scale)
        for key,v in source.items():
            target[key] = target.get(key,0)+scale*v
