"""CONVENTIONS 2.2 sign rules through fastker: reversing k edges multiplies K(lambda;0) by
(-1)^(Dk); a permutation pi of the edge order multiplies it by sgn(pi)^(D+1)."""
import fastker
from check_c3_three_loop import K4SUB, MOM

name, edges, ref = K4SUB[0]
K0 = fastker.kernel(8, edges, MOM, verbose=False)
rev = [(b, a) if i in (0, 3) else (a, b) for i, (a, b) in enumerate(edges)]      # k = 2 -> (+1)
K2 = fastker.kernel(8, rev, MOM, verbose=False)
rev1 = [(b, a) if i == 5 else (a, b) for i, (a, b) in enumerate(edges)]          # k = 1 -> (-1)^3
K1 = fastker.kernel(8, rev1, MOM, verbose=False)
perm = [edges[1], edges[0]] + edges[2:]                                            # transposition -> (-1)^4
Kp = fastker.kernel(8, perm, MOM, verbose=False)
neg = {k: -v for k, v in K0.items()}
print(name, 'K =', {k: str(v) for k, v in K0.items()})
print('reverse 2 edges: same', K2 == K0, '| reverse 1 edge: negated', K1 == neg, '| swap 2 edges (D=3): same', Kp == K0)
