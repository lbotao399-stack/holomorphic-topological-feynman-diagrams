"""Test GL3 x GL2 polynomial functors on one detecting pair of weights."""

from time import monotonic
from diagonal_products import count
from relative import Sector


def balanced(total, size):
    return tuple(total//size + (i<total%size) for i in range(size))


def probe(h, t):
    w = balanced(h,3)+balanced(t,2)
    dims = dict(Sector(w).analyze())[None]
    products = count(w).get(w,0)
    kernel = sum(n-r for n,r in dims.values())
    missing = sum(n for n,r in dims.values())-products
    return w, sum(n for n,r in dims.values()), kernel, missing


if __name__ == "__main__":
    import sys
    total = int(sys.argv[1]) if len(sys.argv)>1 else 8
    for h in range(2,total+1):
        start=monotonic()
        result=probe(h,total-h)
        print(result, round(monotonic()-start,2), flush=True)
