"""Exact bilinear moments after eliminating the last gap on each marked path."""
from array import array
from ctypes import CDLL, c_int, c_void_p
from functools import cache
from math import factorial, lcm, prod
from weakref import finalize
import os

from flint import fmpz_mpoly_ctx


@cache
def native():
    path=os.environ.get('MIXED_SIMPLEX_GRAM_LIB')
    if not path:return None
    lib=CDLL(path)
    make,apply,free=lib.simplex_gram,lib.simplex_gram_row,lib.free_simplex_gram
    make.argtypes=[c_int,c_int,c_int,c_int,c_void_p,c_void_p,c_void_p,c_void_p,
                   c_void_p,c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
    make.restype=c_void_p
    apply.argtypes=[c_void_p,c_int,c_void_p]
    apply.restype=None
    free.argtypes=[c_void_p]
    free.restype=None
    return make,apply,free


@cache
def exact_digits():
    path=os.environ.get('MIXED_SIMPLEX_GRAM_LIB')
    if not path:return None
    library=CDLL(path)
    if not hasattr(library,'digit_simplex_gram'):return None
    f=library.digit_simplex_gram;f.argtypes=[c_void_p,c_int,c_int];f.restype=None
    return f


def pack(rows):
    ends,terms,coefficients=array('i',[0]),array('i'),array('q')
    for row in rows:
        for t,c in row:
            terms.append(t);coefficients.append(c)
        ends.append(len(terms))
    return ends,terms,coefficients


def simplex_variables(dimensions):
    ctx=fmpz_mpoly_ctx.get(tuple('u'+str(k) for k in range(sum(dimensions))))
    variables=ctx.gens()
    substitute=[]
    offset=0
    for size in dimensions:
        group=variables[offset:offset+size]
        substitute.extend((*group,1-sum(group)))
        offset+=size
    return substitute


def sparse_rows(polynomials):
    by_key={}
    polys=[]
    for poly in polynomials:
        terms=[]
        for ex,c in poly.terms():
            key=0
            for power in ex:
                if power>10:raise OverflowError
                key=(key<<5)+int(power)
            terms.append((by_key.setdefault(key,len(by_key)),int(c)))
        polys.append(tuple(terms))
    return tuple(by_key),polys


def gram_table(eh,hh,ew,unique,dimensions):
    if len(eh)**2*len(unique)>12000000 or len(eh)**2*len(ew)>12000000:return None
    ph,pw=pack(hh),pack(unique)
    # Packed keys have five-bit exponent lanes. The three input powers are
    # each at most ten, so no addition can pass a carry to the adjacent lane.
    pairs={i+j for i in eh for j in eh}
    sums=sorted({k+w for k in pairs for w in ew})
    def denominator(key):
        value=1
        for size in reversed(dimensions):
            ex=[]
            for _ in range(size):
                ex.append(key&31);key>>=5
            value*=factorial(size+sum(ex))//prod(map(factorial,ex))
        return value
    denominators=[denominator(k) for k in sums]
    common=lcm(*denominators)
    height=max(sum(abs(c) for k,c in a) for a in hh)
    soft_height=max(sum(abs(c) for k,c in a) for a in unique)
    if (common*height**2*soft_height).bit_length()>126:return None
    lifted=array('Q')
    for deg in denominators:
        value=common//deg
        lifted.extend((value%2**64,value>>64))
    make,apply,free=native()
    arrays=(array('Q',eh),array('Q',ew),*ph,*pw,array('Q',sums),lifted)
    # Monomial keys are sorted only for the sum table. The hard and soft
    # dictionaries retain the column numbering of the sparse polynomial rows.
    storage=make(len(hh),len(eh),len(unique),len(sums),
                 *(x.buffer_info()[0] for x in arrays))
    if not storage:return None
    digit=exact_digits()
    if digit and os.environ.get('MIXED_BINARY_GRAM','1')=='1':
        bits=min(48,52-(max(height,1)**2).bit_length())
        magnitude=(common*soft_height).bit_length()
        if bits>=8 and magnitude<=2*bits:
            digit(storage,bits,max(1,(magnitude+bits-1)//bits))
    return common,storage


def build(hard,weights,dimensions):
    if not native() or not hard or not weights:return None
    substitute=simplex_variables(dimensions)
    try:
        eh,hh=sparse_rows(row[1].compose(*substitute) for row in hard)
        ew,ww=sparse_rows(row[1].compose(*substitute) for row in weights)
        unique=list(dict.fromkeys(ww))
        selected={row:i for i,row in enumerate(unique)}
        cols=[selected[row] for row in ww]
        result=gram_table(eh,hh,ew,unique,dimensions)
    except OverflowError:
        return None
    if result is None:return None
    common,storage=result
    _,apply,free=native()
    raw=array('Q',[0])*(2*len(hh)*len(unique))
    def rows():
        try:
            for i in range(len(hh)):
                apply(storage,i,raw.buffer_info()[0])
                for j in range(i,len(hh)):
                    start=j*2*len(unique)
                    unsigned=(raw[start+2*k]+(raw[start+2*k+1]<<64) for k in range(len(unique)))
                    signed=[v-(1<<128) if v>>127 else v for v in unsigned]
                    yield i,j,[signed[k] for k in cols]
        finally:
            cleanup()
    stream=rows()
    cleanup=finalize(stream,free,storage)
    return common,stream
