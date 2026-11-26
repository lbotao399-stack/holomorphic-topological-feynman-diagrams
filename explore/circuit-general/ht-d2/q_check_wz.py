"""Task 3 checks for Wess-Zumino models in the HT twist (CDG arXiv:2005.00083, ar5iv eqs. (3.32),
(3.34), (3.44)): S = int Psi d' Phi - int W(Phi), {Phi, Psi}_BV = delta, Q Psi = dW/dPhi, Q Phi = 0
on the local operators (holomorphic jets phi^(k) even, psi^(k) odd).

(1) Graph count: an edge needs one antifield terminal (Psi); W vertices have none, so every edge of
    a Q_n graph with only W vertices ends at O: the graph is a star, |I| = n + 1, and the Laman count
    |I| = 2n + 1 forces n = 0.  Checked by enumeration of all tight D = 2 graphs (htk.catalog) with
    every choice of the insertion vertex O, n = 1..4.
(2) Q_0^2 = 0 exactly on all monomials of weight <= wmax (weight = number of letters + derivatives),
    W = Phi^3/3 and W = XYZ.  Hence sum_{k+l=n} Q_k Q_l = 0 for n <= 4 (Q_k = 0 for k >= 1).
Output: q_check_wz.log"""
import itertools
import sys
from fractions import Fraction as Fr
from math import comb

sys.path.insert(0, __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "..", "..", "..", "code"))
from htk import catalog


def graph_count(nmax=4):
    out = {}
    for n in range(1, nmax + 1):
        gs = catalog.tight_graphs(2, n)
        admissible = 0
        for G in gs:
            for O in range(G.nv):
                internal = [e for e in G.edges if O not in e]
                if not internal:
                    admissible += 1
        out[n] = (len(gs), admissible)
    return out


# ---------------------------------------------------------------- graded polynomial algebra
# generator = (field, k): field in fields; parity: 0 (phi-type), 1 (psi-type)

class Alg:
    def __init__(self, phis, psis, Wprime):
        self.phis, self.psis = phis, psis          # names
        self.Wprime = Wprime                        # psi name -> dict monomial (tuple of (phi, k)) -> coeff for W' at order 0

    def parity(self, g):
        return 1 if g[0] in self.psis else 0

    def normal(self, mono):
        """Sort generators (phi's before psi's, then by name, k) with the Koszul sign; psi^2 = 0."""
        mono = list(mono)
        sign = 1
        # bubble sort with signs
        for i in range(len(mono)):
            for j in range(len(mono) - 1 - i):
                a, b = mono[j], mono[j + 1]
                if self.key(a) > self.key(b):
                    mono[j], mono[j + 1] = b, a
                    if self.parity(a) and self.parity(b):
                        sign = -sign
        for j in range(len(mono) - 1):
            if mono[j] == mono[j + 1] and self.parity(mono[j]):
                return None, 0
        return tuple(mono), sign

    def key(self, g):
        return (self.parity(g), g[0], g[1])

    def add(self, P, mono, c):
        m, s = self.normal(mono)
        if m is None or c == 0:
            return
        P[m] = P.get(m, 0) + s * c
        if P[m] == 0:
            del P[m]

    def d_jet(self, poly, times):
        """Apply the holomorphic derivative 'times' times to a polynomial in phi-jets (Leibniz)."""
        cur = dict(poly)
        for _ in range(times):
            nxt = {}
            for mono, c in cur.items():
                for i, g in enumerate(mono):
                    new = list(mono)
                    new[i] = (g[0], g[1] + 1)
                    self.add(nxt, new, c)
            cur = nxt
        return cur

    def Q0_gen(self, g):
        if self.parity(g) == 0:
            return {}
        return self.d_jet(self.Wprime[g[0]], g[1])

    def Q0(self, P):
        """Odd derivation: Q(a_1...a_m) = sum_i (-1)^{|a_1..a_{i-1}|} a_1..Q(a_i)..a_m."""
        out = {}
        for mono, c in P.items():
            sgn = 1
            for i, g in enumerate(mono):
                for m2, c2 in self.Q0_gen(g).items():
                    self.add(out, list(mono[:i]) + list(m2) + list(mono[i + 1:]), sgn * c * c2)
                if self.parity(g):
                    sgn = -sgn
        return out


def monomials(alg, wmax):
    gens = [(f, k) for f in alg.phis + alg.psis for k in range(wmax)]
    out = []
    for size in range(1, wmax + 1):
        for combo in itertools.combinations_with_replacement(gens, size):
            w = sum(1 + g[1] for g in combo)
            if w > wmax:
                continue
            m, s = alg.normal(combo)
            if m is not None:
                out.append(m)
    return sorted(set(out))


def check(alg, wmax):
    bad = 0
    mons = monomials(alg, wmax)
    nonzero_Q = 0
    for m in mons:
        Q1 = alg.Q0({m: 1})
        if Q1:
            nonzero_Q += 1
        Q2 = alg.Q0(Q1)
        if Q2:
            bad += 1
    return len(mons), nonzero_Q, bad


if __name__ == "__main__":
    print("(1) admissible (graph, O) pairs with only W vertices, tight D=2 graphs:")
    for n, (ng, adm) in graph_count(4).items():
        print("   n=%d: %d tight graphs, admissible pairs = %d" % (n, ng, adm))
    one = Alg(["phi"], ["psi"], {"psi": {(("phi", 0), ("phi", 0)): 1}})          # W = phi^3/3
    xyz = Alg(["X", "Y", "Z"], ["pX", "pY", "pZ"], {"pX": {(("Y", 0), ("Z", 0)): 1},
                                                    "pY": {(("X", 0), ("Z", 0)): 1},
                                                    "pZ": {(("X", 0), ("Y", 0)): 1}})
    for name, alg, w in (("W = Phi^3/3", one, 6), ("W = XYZ", xyz, 4)):
        nm, nq, bad = check(alg, w)
        print("(2) %s: %d monomials of weight <= %d, Q_0 nonzero on %d, Q_0^2 != 0 on %d" % (name, nm, w, nq, bad))
    # example
    m = (("phi", 1), ("psi", 0), ("psi", 2))
    print("    example Q_0(dphi psi d^2psi) =", one.Q0({one.normal(m)[0]: one.normal(m)[1]}))
