"""Optional native contraction of the marked-position factorial moments."""
from array import array
from ctypes import CDLL, c_int, c_void_p
from functools import cache
from math import factorial, prod
from weakref import finalize
import os


@cache
def native_kernel():
    path = os.environ.get('MIXED_PRISM_MOMENTS_LIB')
    if not path:
        return None
    function = CDLL(path).ordered_moments
    function.argtypes = [c_int,c_int,c_void_p,c_void_p,c_int,c_void_p,c_void_p,c_void_p,c_void_p]
    function.restype = c_int
    return function


@cache
def native_table():
    path=os.environ.get('MIXED_PRISM_MOMENTS_LIB')
    if not path:return None
    library=CDLL(path)
    if not hasattr(library,'ordered_moment_table'):return None
    make,apply,free=library.ordered_moment_table,library.table_moments,library.free_ordered_moment_table
    make.argtypes=[c_int,c_int,c_void_p,c_void_p,c_void_p]
    make.restype=c_void_p
    apply.argtypes=[c_void_p,c_int,c_void_p,c_void_p,c_void_p]
    apply.restype=c_int
    free.argtypes=[c_void_p]
    free.restype=None
    return make,apply,free


@cache
def native_product():
    path=os.environ.get('MIXED_PRISM_MOMENTS_LIB')
    if not path:return None
    library=CDLL(path)
    if not hasattr(library,'product_moments'):return None
    apply=library.product_moments
    apply.argtypes=[c_void_p,c_int,c_void_p,c_void_p,c_int,c_void_p,c_void_p,c_void_p]
    apply.restype=c_int
    return apply


def signed_values(out):
    values=[lo+(hi<<64) for lo,hi in zip(out[::2],out[1::2])]
    return [v-(1<<128) if v>>127 else v for v in values]


def packed(polynomials):
    powers,coefficients,ends = array('B'),array('q'),array('i',[0])
    for p in polynomials:
        for ex,c in p.terms():
            powers.extend(map(int,ex))
            coefficients.append(int(c))
        ends.append(len(coefficients))
    return powers,coefficients,ends


class PositionMoments:
    def __init__(self, weights):
        self.weights = weights
        self.kernel = native_kernel()
        self.table=None
        if self.kernel:
            try:
                self.powers,self.coefficients,self.ends = packed(weights)
                functions=native_table()
                if functions and self.coefficients:
                    make,self.apply,free=functions
                    self.table=make(weights[0].context().nvars(),len(weights),self.ends.buffer_info()[0],
                                    self.powers.buffer_info()[0],self.coefficients.buffer_info()[0])
                    if self.table:finalize(self,free,self.table)
            except OverflowError:
                self.kernel = None

    def factor(self,poly):
        if not self.table or not native_product():return None
        powers,coeffs=array('Q'),array('q')
        maxima=[0]*poly.context().nvars()
        try:
            for ex,c in poly.terms():
                key=0
                for i,degree in enumerate(ex):
                    degree=int(degree)
                    if degree>20:return None
                    key=(key<<5)|degree
                    maxima[i]=max(maxima[i],degree)
                powers.append(key)
                coeffs.append(int(c))
        except OverflowError:
            return None
        return powers,coeffs,maxima

    def product(self,x,y,left,right):
        if left is not None and right is not None:
            px,ax,dx=left
            py,ay,dy=right
            if all(a+b<=20 for a,b in zip(dx,dy)):
                out=array('Q',[0])*(2*len(self.weights))
                if native_product()(self.table,len(ax),px.buffer_info()[0],ax.buffer_info()[0],
                                    len(ay),py.buffer_info()[0],ay.buffer_info()[0],out.buffer_info()[0]):
                    return signed_values(out)
        return self(x*y)

    def __call__(self, local):
        if self.kernel:
            try:
                powers,coefficients,_ = packed([local])
                n = local.context().nvars()
                out = array('Q',[0])*(2*len(self.weights))
                if self.table:
                    ok=self.apply(self.table,len(coefficients),powers.buffer_info()[0],coefficients.buffer_info()[0],
                                  out.buffer_info()[0])
                else:
                    ok=self.kernel(n,len(coefficients),powers.buffer_info()[0],coefficients.buffer_info()[0],
                                   len(self.weights),self.ends.buffer_info()[0],self.powers.buffer_info()[0],
                                   self.coefficients.buffer_info()[0],out.buffer_info()[0])
                if ok:
                    return signed_values(out)
            except OverflowError:
                pass
        return [sum(int(c)*prod(factorial(int(j)) for j in ex) for ex,c in (local*q).terms())
                for q in self.weights]
