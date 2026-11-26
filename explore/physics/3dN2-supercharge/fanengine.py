"""Q_n on undifferentiated operators (z = 0), cubic theories (kappa = 0), from canonical graphs:
  fan(n):            htk graph.fan(n), path 0..n, root R = n+1, K = s_n C_n (l_0 l_n)^{n/2}
  tadpole(n+1, 0):   htk tadpole_graph(n+1, 0), cycle 0..n, q = 0, spokes (v, R), K = t_n l_0^n
Q_n^{shape} O = (1/|Aut|) sum over patterns on the canonical graph (vertex i of the graph = i-th copy
of I), |Aut| = 2 for both shapes (path reversal; reflection of the cycle fixing q).
Kernel normalization: values divided by pi^{n/2}."""
from fractions import Fraction as Fr
from itertools import product, permutations
from ht3d import deriv_dist, padd


def fan_graph(n):
    R = n + 1
    edges = []
    for i in range(n):
        edges += [(i, i + 1), (i, R)]
    edges.append((n, R))
    return edges, R


def tadpole_graph(c):
    n = c - 1
    R = n + 1
    edges = [(i, (i + 1) % c) for i in range(c)]
    edges += [(v, R) for v in range(1, n + 1)]
    edges = sorted((min(a, b), max(a, b)) for a, b in edges)
    return edges, R


class ShapeEngine:
    def __init__(self, th, n, shape, kval):
        """kval(exps) -> Fraction: kernel coefficient of prod_v l_v^{exps[v]} (z = 0, divided by pi^{n/2})"""
        self.th, self.A, self.n = th, th.A, n
        if shape == 'fan':
            self.edges, self.R = fan_graph(n)
        else:
            self.edges, self.R = tadpole_graph(n + 1)
        self.kval = kval
        self.cache = {}

    def mono(self, O):
        if O in self.cache:
            return self.cache[O]
        th, A, n, R = self.th, self.A, self.n, self.R
        edges = self.edges
        Il = th.I_list
        out = {}
        nv = n + 1
        # neighbours of each interaction vertex among smaller vertices, and spoke flags
        spoke = [any(e == (v, R) for e in edges) for v in range(nv)]
        lower = [[u for (u, w) in edges if w == v] for v in range(nv)]

        def compat(m1, m2):
            return any(th.eta(x, y) for x in m1 for y in m2)

        def rec(v, chosen):
            if v == nv:
                yield list(chosen)
                return
            for (m, c) in Il:
                if spoke[v] and not compat(m, O):
                    continue
                if any(not compat(chosen[u][0], m) for u in lower[v]):
                    continue
                chosen.append((m, c))
                yield from rec(v + 1, chosen)
                chosen.pop()

        for chosen in rec(0, []):
            mons = [m for m, _ in chosen] + [O]
            coeffs = [c for _, c in chosen]
            cands = []
            ok = True
            for (v, w) in edges:
                cl = [(i, j) for i, x in enumerate(mons[v]) for j, y in enumerate(mons[w]) if th.eta(x, y)]
                if not cl:
                    ok = False
                    break
                cands.append(cl)
            if not ok:
                continue
            for choice in product(*cands):
                used = set()
                good = True
                for (v, w), (i, j) in zip(edges, choice):
                    if (v, i) in used or (w, j) in used:
                        good = False
                        break
                    used.add((v, i)); used.add((w, j))
                if not good:
                    continue
                # Wick sign
                letters, offs = [], []
                for m in mons:
                    offs.append(len(letters)); letters += list(m)
                order, eta = [], 1
                gused = set()
                for (v, w), (i, j) in zip(edges, choice):
                    gi, gj = offs[v] + i, offs[w] + j
                    order += [gi, gj]; gused.add(gi); gused.add(gj)
                    eta *= th.eta(letters[gi], letters[gj])
                surv = [g for g in range(len(letters)) if g not in gused]
                sign = A.perm_sign(letters, order + surv)
                base = Fr(sign * eta)
                for c in coeffs:
                    base *= c
                groups = [[letters[offs[v] + i] for i in range(len(mons[v])) if offs[v] + i not in gused]
                          for v in range(len(mons))]
                # kernel: only monomials whose lambda-degree at v matches survivors (any distribution)
                for exps, kc in self.kval_terms(groups):
                    parts = []
                    okp = True
                    for v in range(nv):
                        dd = deriv_dist(groups[v], exps[v])
                        if not dd:
                            okp = False
                            break
                        parts.append(dd)
                    if not okp:
                        continue
                    for combo in product(*parts):
                        mult = 1
                        lets = []
                        for mlt, ls in combo:
                            mult *= mlt; lets += list(ls)
                        lets += groups[nv]
                        s, mm = A.canon(lets)
                        if s:
                            padd(out, mm, base * kc * mult * s)
        out = {m: c / 2 for m, c in out.items() if c}
        self.cache[O] = out
        return out

    def kval_terms(self, groups):
        return self.kval(groups)

    def __call__(self, P):
        out = {}
        for O, c in P.items():
            for m, v in self.mono(O).items():
                padd(out, m, c * v)
        return out


def fan_kval(n, s_n_Cn):
    """K = s_n C_n (l_0 l_n)^{n/2}: exps with l_0^{n/2} l_n^{n/2}"""
    h = n // 2

    def kv(groups):
        e = [0] * (n + 1)
        e[0] = h; e[n] = h
        return [(tuple(e), s_n_Cn)]
    return kv


def tad_kval(n, t_n):
    def kv(groups):
        e = [0] * (n + 1)
        e[0] = n
        return [(tuple(e), t_n)]
    return kv
