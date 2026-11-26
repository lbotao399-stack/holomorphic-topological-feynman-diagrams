"""(a) SL2-invariant polynomials of degree 4 and bidegree (2,2) in k free vectors of C^2 (k = 2, 3, 4):
    dim ker E, E = sum_i q_i^1 d/dq_i^2, on bidegree-(2,2) polynomials (= SL2 invariants: weight-0 vectors killed by
    the raising operator of a finite-dimensional sl2-module).  Candidates: products of two brackets; rank printed.
(b) [q_s0 q_s1] = sgn(s) [q0 q1] modulo q0 + q1 + q2 = 0 for all s in S3.
(c) degree 2 (bidegree (1,1)) in one free vector: dim ker E (SL2-invariant 2-cochain polynomials of J_1 degree)."""
import itertools
import sympy as sp

def run(k, deg=(2, 2)):
    q = [(sp.Symbol('a%d' % i), sp.Symbol('b%d' % i)) for i in range(k)]
    first = [x[0] for x in q]
    second = [x[1] for x in q]
    m1 = [sp.Mul(*[first[i] ** e[i] for i in range(k)]) for e in itertools.product(range(deg[0] + 1), repeat=k) if sum(e) == deg[0]]
    m2 = [sp.Mul(*[second[i] ** e[i] for i in range(k)]) for e in itertools.product(range(deg[1] + 1), repeat=k) if sum(e) == deg[1]]
    monos = [a * b for a in m1 for b in m2]
    E = lambda P: sp.expand(sum(q[i][0] * sp.diff(P, q[i][1]) for i in range(k)))
    gens = first + second
    images = [sp.Poly(E(m), *gens) for m in monos]
    keys = sorted({mm for P in images for mm in P.as_dict()})
    M = sp.Matrix([[P.as_dict().get(kk, 0) for P in images] for kk in keys]) if keys else sp.zeros(1, len(monos))
    dker = len(monos) - M.rank()
    # candidates
    br = lambda i, j: q[i][0] * q[j][1] - q[i][1] * q[j][0]
    pairs = list(itertools.combinations(range(k), 2))
    if sum(deg) == 4:
        cands = [sp.expand(br(*pairs[a]) * br(*pairs[b])) for a, b in itertools.combinations_with_replacement(range(len(pairs)), 2)]
    elif sum(deg) == 2:
        cands = [sp.expand(br(*p)) for p in pairs]
    else:
        cands = []
    if cands:
        cp = [sp.Poly(c, *gens) for c in cands]
        ck = sorted({mm for P in cp for mm in P.as_dict()})
        rc = sp.Matrix([[P.as_dict().get(kk, 0) for P in cp] for kk in ck]).rank()
        killed = all(E(c) == 0 for c in cands)
    else:
        rc, killed = 0, True
    return len(monos), dker, len(cands), rc, killed

for k in (2, 3, 4):
    nm, dk, nc, rc, kl = run(k)
    print("free vectors %d, bidegree (2,2): monomials %d, dim ker E = %d ; bracket products %d, rank %d, killed by E: %s"
          % (k, nm, dk, nc, rc, kl))
nm, dk, nc, rc, kl = run(1, (1, 1))
print("free vectors 1, bidegree (1,1): monomials %d, dim ker E = %d" % (nm, dk))
nm, dk, nc, rc, kl = run(2, (1, 1))
print("free vectors 2, bidegree (1,1): monomials %d, dim ker E = %d ; brackets %d, rank %d, killed: %s" % (nm, dk, nc, rc, kl))
# (b)
a = [sp.symbols('x%d y%d' % (i, i)) for i in range(2)]
qs = [a[0], a[1], (-a[0][0] - a[1][0], -a[0][1] - a[1][1])]
br = lambda p, r: p[0] * r[1] - p[1] * r[0]
ok = True
for s in itertools.permutations(range(3)):
    sg = 1
    for i in range(3):
        for j in range(i + 1, 3):
            if s[i] > s[j]:
                sg = -sg
    ok &= sp.expand(br(qs[s[0]], qs[s[1]]) - sg * br(qs[0], qs[1])) == 0
print("[q_s0 q_s1] = sgn(s)[q0 q1] mod conservation, all s in S3:", ok)
