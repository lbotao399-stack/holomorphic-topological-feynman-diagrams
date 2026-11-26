"""Contract the five exterior path moments without constructing the full exterior polynomial."""
import ctypes
from collections import defaultdict
from fractions import Fraction
from functools import cache,lru_cache
from hashlib import sha256
from itertools import combinations_with_replacement
from math import factorial,gcd,lcm
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from flint import fmpz_mpoly_ctx
from quartic_seven_plane_single import cuts
from quartic_seven_plane_single_projection import Bubble,Q,SCALE,orbital,split_variables


@cache
def native_fold():
    source = Path(__file__).with_suffix('.cpp')
    target = Path(tempfile.gettempdir())/('mixed-bare-kite-'+sha256(source.read_bytes()).hexdigest()[:16]+'.so')
    if not target.exists():
        fd,name = tempfile.mkstemp(dir=target.parent,suffix='.so')
        os.close(fd)
        try:
            subprocess.run(['c++','-O3','-std=c++17','-shared','-fPIC',str(source),'-o',name],check=True)
            os.replace(name,target)
        finally:
            Path(name).unlink(missing_ok=True)
    fn = ctypes.CDLL(str(target)).fold
    p,u = ctypes.c_void_p,ctypes.c_uint
    fn.argtypes = [u,p,p,u,p,p,p,u,p,p,p]
    fn.restype = ctypes.c_int
    return fn


def packed(poly):
    return (np.fromiter((int(ex[0]) for ex in poly.monoms()),np.uint64,len(poly)),
            np.fromiter((int(c) for c in poly.coeffs()),np.int64,len(poly)))


def sparse_split(poly,nleft):
    """Exact sparse row elimination, packing exponent vectors into base sixteen."""
    from quartic_seven_plane_double import primitive

    nright = poly.context().nvars()-nleft
    slices = defaultdict(lambda:defaultdict(dict))
    for mon,c in poly.terms():
        ex = tuple(map(int,mon))
        assert max(ex)<16
        left = sum(k<<(4*i) for i,k in enumerate(ex[:nleft]))
        right = sum(k<<(4*i) for i,k in enumerate(ex[nleft:]))
        slices[sum(ex[:nleft])][right][left] = int(c)

    def unpack(key,size):
        return tuple((key>>(4*i))&15 for i in range(size))

    for degree in sorted(slices):
        entries = slices.pop(degree)
        basis = []
        for right in sorted(entries,key=lambda r:(len(entries[r]),r)):
            row = entries.pop(right)
            denominator = 1
            for pivot,left,column in basis:
                value = row.get(pivot,0)
                if not value:
                    continue
                leading = left[pivot]
                column[right] = Fraction(value,leading*denominator)
                common = gcd(value,leading)
                factor,multiple = leading//common,value//common
                if factor!=1:
                    row = {k:v*factor for k,v in row.items()}
                    denominator *= factor
                for k,v in left.items():
                    new = row.get(k,0)-multiple*v
                    if new:
                        row[k] = new
                    else:
                        row.pop(k,None)
            if row:
                pivot = min(row,key=lambda k:(abs(row[k]),k))
                content = gcd(*row.values())*(1 if row[pivot]>0 else -1)
                basis.append((pivot,{k:v//content for k,v in row.items()},
                              {right:Fraction(content,denominator)}))
        for _,left,column in basis:
            l = tuple(sorted((unpack(k,nleft),v) for k,v in left.items()))
            r,factor = primitive(sorted((unpack(k,nright),v) for k,v in column.items()))
            yield l,r,degree,factor


def project(data,squares,voltage,check_star=False,statistics=None,sparse=False):
    paths,split = data[5:]
    nleft = sum(map(len,paths[:4]))
    sizes = np.array([len(path) for path in paths[4:]],dtype=np.uint32)
    right = fmpz_mpoly_ctx.get(('packed',))
    bubble = Bubble(nleft,split)
    polynomials,ids = [],{}
    def parts(poly,max_degree):
        result = []
        separate = sparse_split if sparse else split_variables
        for l,r,degree,factor in separate(poly,nleft):
            if r not in ids:
                encoded = {}
                for ex,value in r:
                    if max(ex)>max_degree:
                        raise ArithmeticError('Exterior segment exceeds its encoding degree bound')
                    encoded[(sum(n<<(4*i) for i,n in enumerate(ex)),)] = value
                ids[r] = len(polynomials)
                polynomials.append(right.from_dict(encoded))
            result.append((bubble.register(l),ids[r],degree,factor))
        return result
    weights = defaultdict(Fraction)
    singular = defaultdict(Fraction)
    for root,position in squares:
        # Squaring two degree-three path factors and adding the voltage and position
        # powers stays below sixteen in every digit of the sparse Kronecker product.
        pieces = parts(root,3)
        third = parts(position*voltage,5)
        for i,j in combinations_with_replacement(range(len(pieces)),2):
            a,r,g,c = pieces[i]
            b,s,h,d = pieces[j]
            pair = tuple(sorted((r,s)))
            for k,t,shift,e in third:
                local,lost = bubble.project(a,b,k)
                coefficient = c*d*e*(1+(i!=j))
                if local:
                    weights[g+h+shift+nleft-Q,*pair,t] += coefficient*local
                for index,value in lost:
                    singular[index,*pair,t] += coefficient*value
    # Only the combined coefficients of locally nonintegrable bubble monomials can enter.
    if any(singular.values()):
        from flint import fmpq,fmpq_mpoly_ctx
        qq = fmpq_mpoly_ctx.get(right.names())
        bad = defaultdict(lambda:qq.constant(0))
        for (index,r,s,t),w in singular.items():
            if w:
                poly = polynomials[r]*polynomials[s]*polynomials[t]
                bad[index] += qq.from_dict(poly.to_dict())*fmpq(w.numerator,w.denominator)
        assert not any(bad.values())

    @lru_cache(maxsize=6)
    def square_pair(r,s):
        return packed(polynomials[r]*polynomials[s])
    @cache
    def third_factor(t):
        return packed(polynomials[t])
    capacity = 32768
    keys = np.empty(capacity,dtype=np.uint32)
    low = np.empty(capacity,dtype=np.uint64)
    high = np.empty(capacity,dtype=np.int64)
    # Each external power sum is integral after one common denominator is
    # cleared. The very many exterior terms then require only integer adds.
    denominator = lcm(*(w.denominator for w in weights.values()))
    folded = defaultdict(int)
    singles = [i for i,n in enumerate(sizes) if n==1]
    exterior_products = 0
    exterior_rows = 0
    for (power,r,s,t),weight in sorted(weights.items(),key=lambda item:item[0][1:]):
        if not weight: continue
        multiplier = weight.numerator*(denominator//weight.denominator)
        lm,lc = square_pair(r,s)
        rm,rc = third_factor(t)
        n = native_fold()(len(lm),lm.ctypes.data,lc.ctypes.data,len(rm),rm.ctypes.data,rc.ctypes.data,
                          sizes.ctypes.data,capacity,keys.ctypes.data,low.ctypes.data,high.ctypes.data)
        if n<0: raise ArithmeticError(f'Exterior integer contraction exceeded its exact bounds ({n})')
        exterior_products += 1
        exterior_rows += n
        if n:
            key = int(keys[0])
            assert power+sum((key>>(4*j))&15 for j in range(5))+sum(map(int,sizes))==3*Q
        for i in range(n):
            integer = int(low[i])+(int(high[i])<<64)
            if integer: folded[int(keys[i])] += integer*multiplier
    normalized = {}
    for key,weight in folded.items():
        if not weight: continue
        degrees = [(key>>(4*j))&15 for j in range(5)]
        rest = tuple(d+int(k) for d,k in zip(degrees,sizes))
        power = 3*Q-sum(rest)
        factor = factorial(power-1)
        for j in singles: factor *= factorial(degrees[j])
        normalized[(power,*rest)] = Fraction(weight*factor,denominator*SCALE*factorial(Q-1))
    reduced = defaultdict(Fraction)
    for signature,weight in normalized.items():
        if weight: reduced[orbital(signature)] += weight
    value = cuts({key:w for key,w in reduced.items() if w})
    if statistics is not None:
        statistics.update(bubble_triples=len(bubble.triples),exterior_products=exterior_products,
                          exterior_rows=exterior_rows)
    if check_star:
        from quartic_seven_plane_single_star import star_integral
        assert value==(star_integral({key:w for key,w in normalized.items() if w}),0)
    return value
