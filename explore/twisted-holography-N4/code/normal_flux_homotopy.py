"""Smooth polynomial flux modes and a homotopy commuting with holomorphic derivatives."""
from functools import lru_cache
from itertools import combinations_with_replacement, permutations
from math import factorial

import sympy as sp


W = sp.symbols("w1:4")
R = sp.symbols("r1:4")
G = sum(w*r for w, r in zip(W, R))


def trace(p):
    return sum(sp.diff(p, w, r) for w, r in zip(W, R)).expand()


@lru_cache(None)
def harmonics(p, k, m):
    """Decompose a bihomogeneous polynomial as a sum of powers of G times trace-free terms."""
    top = derivative = p
    for j in range(1, min(k, m)+1):
        derivative = trace(derivative)
        if not derivative:
            break
        scale = sp.Rational((-1)**j*factorial(k+m+1-j), factorial(j)*factorial(k+m+1))
        top += scale*G**j*derivative
    top = top.expand()
    assert trace(top) == 0
    terms = [(top, k, m)]
    if top != p:
        quotient = sp.div(p-top, G, *W, *R)[0]
        assert (p-top-G*quotient).expand() == 0
        terms += harmonics(quotient, k-1, m-1)
    return terms


def split(p):
    groups = {}
    for powers, value in sp.Poly(p, *W, *R).terms():
        if value:
            degree = sum(powers[:3]), sum(powers[3:])
            monomial = sp.prod(x**n for x, n in zip((*W, *R), powers))
            groups[degree] = groups.get(degree, 0) + value*monomial
    return groups


def transfer(p):
    """Return coefficients of b_i in the primitive and of U in the Laurent projection."""
    primitive = sp.zeros(3, 1)
    harmonic = sp.S.Zero
    for (k, m), group in split(p).items():
        for term, a, b in harmonics(group.expand(), k, m):
            if a:
                primitive += sp.Rational(2, a*(b+2))*sp.Matrix([term.diff(w) for w in W])
            else:
                harmonic += term
    return primitive, harmonic


def transfer_one(form):
    curl = sp.Matrix([sp.diff(form[(i+2)%3], W[(i+1)%3]) -
                      sp.diff(form[(i+1)%3], W[(i+2)%3]) for i in range(3)])
    scalar = sp.S.Zero
    for (k, m), group in split(sp.Matrix(R).dot(curl)).items():
        for term, a, b in harmonics(group.expand(), k, m):
            if b:
                scalar -= term / (b*(a+2))
    return scalar.expand()


def scalar_dolbeault(p):
    return sp.Matrix(W).cross(sp.Matrix([sp.diff(p, r) for r in R]))


def linking_residue(p):
    result = sp.S.Zero
    for powers, coefficient in sp.Poly(p, *W, *R).terms():
        if powers[:3] == powers[3:]:
            degree = sum(powers[:3])
            result += coefficient*2*sp.prod(factorial(n) for n in powers[:3])/factorial(degree+2)
    return result


def wedge_one(first, second):
    return sp.Matrix(R).dot(first.cross(second))/2


def holomorphic(p, i, power):
    """Coefficient derivative on U (power three) or b_j (power two)."""
    return (sp.diff(p, W[i]) - R[i]*(power*p + sum(r*sp.diff(p, r) for r in R))).expand()


def dolbeault(primitive):
    out = sp.S.Zero
    for i, component in enumerate(primitive):
        for (_, m), p in split(component).items():
            out += ((m+2)*W[i]*p-sp.diff(p, R[i]))/2
    return out.expand()


def on_surface(p):
    return sp.div(sp.expand(p), G-1, *W, *R)[1] == 0


def regular(p):
    result = sp.S.Zero
    for (k, m), group in split(p).items():
        result += sum(term for term, a, b in harmonics(group.expand(), k, m) if b == 0)
    return result


def alternated_trace_operator(p):
    result = sp.S.Zero
    for i, j, k in permutations(range(3)):
        first = transfer(W[k]*p)[0] - W[k]*transfer(p)[0]
        second = transfer_one(W[j]*first) - W[j]*transfer_one(first)
        sign = 1 if (j-i)%3 == 1 else -1
        result += sign*(regular(W[i]*second) - W[i]*regular(second))
    return result.expand()


def check():
    fixtures = [w*r for w in W for r in R]
    fixtures += [W[0]*W[1]*R[1]*R[2], W[1]**2*R[1]**3, G**2*W[0]*R[2],
                 W[0]**3, R[0]*R[1]**2, 1+W[2]*R[2]]
    one_forms = [sp.Matrix([1, W[0]*R[1], W[1]*R[0]**2]),
                 sp.Matrix([W[2]*W[0]*R[0], W[1]*W[0]*R[2], W[0]+R[1]]),
                 sp.Matrix(R)*(W[0]+R[1])]
    for p in fixtures:
        h, projected = transfer(p)
        assert on_surface(dolbeault(h)+projected-p)
        assert on_surface(transfer_one(h))
        for form in one_forms:
            assert linking_residue(wedge_one(h, form)) == linking_residue(p*transfer_one(form))
        for i in range(3):
            dh, dp = transfer(holomorphic(p, i, 3))
            difference = dh-h.applyfunc(lambda q: holomorphic(q, i, 2))
            assert all(on_surface(q) for q in sp.Matrix(R).cross(difference))
            assert on_surface(dp-holomorphic(projected, i, 3))
    for form in one_forms:
        scalar = transfer_one(form)
        upper, projected = transfer(dolbeault(form))
        assert on_surface(projected)
        residue = sp.Matrix(R).cross(form - scalar_dolbeault(scalar) - upper)
        assert all(on_surface(p) for p in residue)
        for i in range(3):
            assert on_surface(transfer_one(form.applyfunc(lambda p: holomorphic(p, i, 2)))
                              - holomorphic(scalar, i, 0))
    for p in [1, W[0], W[0]*W[1], W[2]**3, R[1], W[1]*R[1], W[1]*W[2]*R[1], G]:
        assert on_surface(alternated_trace_operator(p) + linking_residue(p))
    # Pure holomorphic derivatives of the three seed pairs have one radial relation.
    for m in range(6):
        words = list(combinations_with_replacement(range(3), m))
        higher = list(combinations_with_replacement(range(3), m+1))
        rows = {(j, word): j*len(higher)+i for j in range(3) for i, word in enumerate(higher)}
        matrix = sp.MutableSparseMatrix(3*len(higher), 3*len(words), {})
        for j in range(3):
            for i, word in enumerate(words):
                for sign, step in ((1, 1), (-1, 2)):
                    root, component = (j+step)%3, (j-step)%3
                    row = rows[component, tuple(sorted((*word, root)))]
                    matrix[row, j*len(words)+i] = sign
        assert matrix.rank() == (m+1)*(m+3)
    print("Smooth flux contraction and Laurent projection commute with all three coefficient derivatives.")
    print("The alternating three coordinate insertions with unit flux have supertrace -1.")


if __name__ == '__main__':
    check()
