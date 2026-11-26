"""Lemma (closedness of Delta^n T): for five momenta p_0..p_4 in C^3 with sum 0 put
  c^(n)_{ij} = (-1)^{i+j} det(p_k,p_l,p_m)^n,   {k<l<m} = complement of {i<j},  c_{ji} = -c_{ij}.
d(Delta^n T)(y_0..y_4) = (-1)^n sum_{i<j} c^(n)_{ij} T([x_i,x_j],x_k,x_l,x_m)  (T symmetric, n odd; T alternating, n even: see ce_check.py).
For invariant T the vectors T_{ij} = T([x_i,x_j],...) are antisymmetric with zero row sums (cycle space Z of K_5, dim 6).
Hence d(Delta^n T) = <c^(n), T> vanishes identically for every invariant T iff c^(n) is a gradient (c_ij = phi_i - phi_j).
Exact symbolic check (sympy) of: c^(1) gradient; c^(n), n = 2..6, not a gradient (nonzero triangle sum)."""
import sympy as sp
from itertools import combinations
P = [[sp.Symbol('p%d_%d' % (v, a)) for a in range(3)] for v in range(4)]
P.append([-sum(P[v][a] for v in range(4)) for a in range(3)])
def det(a, b, c):
    return sp.Matrix([a, b, c]).det()
delta = {}
for i, j in combinations(range(5), 2):
    k, l, m = [v for v in range(5) if v not in (i, j)]
    delta[(i, j)] = sp.expand(det(P[k], P[l], P[m]))
def c(n, i, j):
    if i < j:
        return (-1)**(i + j) * delta[(i, j)]**n
    return -(-1)**(i + j) * delta[(j, i)]**n
for n in range(1, 7):
    tri = [sp.expand(c(n, 0, 1) + c(n, 1, 2) + c(n, 2, 0)), sp.expand(c(n, 0, 1) + c(n, 1, 3) + c(n, 3, 0)),
           sp.expand(c(n, 1, 2) + c(n, 2, 4) + c(n, 4, 1))]
    allzero = all(sp.expand(c(n, a, b) + c(n, b, cc) + c(n, cc, a)) == 0 for a, b, cc in combinations(range(5), 3))
    print('n=%d: all 10 triangle sums c_ab + c_bc + c_ca vanish: %s' % (n, allzero))
    if n == 1:
        ok = all(sp.expand(c(1, i, j) - (c(1, 0, j) - c(1, 0, i))) == 0 for i, j in combinations(range(1, 5), 2))
        print('   n=1: c_ij = c_0j - c_0i for all 1<=i<j<=4:', ok)
    else:
        print('   n=%d: triangle (0,1,2) sum nonzero: %s' % (n, tri[0] != 0))

# even n: alternating T, T_ij := (-1)^{i+j} T([x_i,x_j], complement sorted) is symmetric with zero row sums;
# d(Delta^n T) = +- sum_{i<j} delta_ij^n T_ij; closed for all invariant alternating T iff delta^n_ij = phi_i + phi_j.
# Test of the phi_i + phi_j form: e_01 + e_23 - e_02 - e_13 = 0 and e_01 + e_34 - e_03 - e_14 = 0 are necessary.
print()
for n in (0, 2, 4, 6):
    e = {k: v**n for k, v in delta.items()}
    t1 = sp.expand(e[(0, 1)] + e[(2, 3)] - e[(0, 2)] - e[(1, 3)])
    # sufficiency at n=0 is trivial (phi = 1/2)
    print('even n=%d: e_01+e_23-e_02-e_13 == 0: %s' % (n, t1 == 0))
