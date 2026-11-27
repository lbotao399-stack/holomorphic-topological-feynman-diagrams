"""Kernel polynomials I_Gamma(lambda) (no input jets) of the low-order rooted graphs, rgraphs labelling.
I = (-1)^{n(n+1)/2} K_Gamma (htk), lambda_v = output momentum at vertex v.  Log: kernels_closed.log
usage: python3 kernels_closed.py [H n] ..."""
import sys
from fractions import Fraction as Fr

import sympy as sp

import rgraphs
import kern_h
from kern_h import GraphKernel, P


def poly_of(gk):
    lexp, coef = gk.kernel_symbolic(((),) * len(gk.spokes))
    syms = {nm: sp.Symbol(nm) for nm in gk.lam_names}
    expr = 0
    for t in range(lexp.shape[0]):
        c = kern_h.ratrec(int(coef[t]))
        if c is None:
            return None
        term = sp.Rational(c.numerator, c.denominator)
        for i, nm in enumerate(gk.lam_names):
            term *= syms[nm] ** int(lexp[t, i])
        expr += term
    return expr


def main(pairs):
    for H, n in pairs:
        gs = rgraphs.rooted_trees(H, n)
        uni = rgraphs.rooted_unicyclic(H, n) if H >= 3 else []
        for k, g in enumerate(gs + uni, 1):
            gk = GraphKernel(H, g['nb'], g['edges'])
            e = poly_of(gk)
            fe = sp.factor(e) if e is not None else 'rational reconstruction failed'
            nterms = len(sp.Add.make_args(sp.expand(e))) if e is not None else '-'
            s = str(fe)
            if len(s) > 400:
                s = s[:400] + ' ...'
            print('H=%d n=%d #%d %s %s |Aut|=%d outputs %s: I = %s  [%s monomials]' % (
                H, n, k, g['kind'], g['edges'], g['aut'], gk.outputs, s, nterms), flush=True)


if __name__ == '__main__':
    a = [int(x) for x in sys.argv[1:]]
    pairs = list(zip(a[::2], a[1::2])) or [(1, 1), (2, 1), (2, 2), (2, 3), (3, 1), (3, 2), (4, 1), (4, 2), (5, 1)]
    main(pairs)
