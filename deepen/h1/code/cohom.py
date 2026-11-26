"""Cohomology of V = M / sum_u D_u M (translation-invariant local functionals, jet model) at fixed
weight w = (w_tau, w_zb1, w_zb2), holomorphic bidegree (a, b), ghost number g; and of the holomorphic
CE complex Q = C*_red(g[[z1,z2]]) / (D_z1, D_z2) at bidegree (a, b), CE degree n (holo=True).

dim H^g = dim V^g - rk Qbar^g - rk Qbar^{g-1},
dim V^g = dim M^g - rk D^g,     rk Qbar^g = rk[Q M^g | D M^{g+1}] - rk D^{g+1}.
Every rank: mod P1, mod P2 and exact over Q (fraction-free elimination); the three are printed.
"""
import sys, time
sys.path.insert(0, __file__.rsplit("/", 1)[0])
import jets as J
import sla


def budget_of(w, ab):
    return (w[0] + 1, w[1] + 1, w[2] + 1, ab[0], ab[1])


class Piece:
    def __init__(self, th, w, ab, holo=False, charge=None):
        self.th, self.w, self.ab, self.holo, self.charge = th, w, ab, holo, charge
        self.B = (0, 0, 0, ab[0], ab[1]) if holo else budget_of(w, ab)
        cf = getattr(th, "charge_fn", J.sl2_charge) if charge is not None else None
        self.cf = cf
        self.M = J.enumerate_monos(th.dim, self.B, holo_only=holo, charge=charge, charge_fn=cf)
        # lower budgets for total derivatives
        self.lower = []
        dirs = [(3, 1), (4, 2)] if holo else [(0, 0), (1, 3), (2, 4), (3, 1), (4, 2)]
        # (index in budget, direction u)
        for bi, u in dirs:
            if self.B[bi] == 0:
                continue
            Bl = list(self.B); Bl[bi] -= 1
            Ml = J.enumerate_monos(th.dim, tuple(Bl), holo_only=holo, charge=charge, charge_fn=cf)
            self.lower.append((u, Ml))
        self._DM = {}

    def Mg(self, g):
        return self.M.get(g, [])

    def DM(self, g):
        if g not in self._DM:
            vecs = []
            for u, Ml in self.lower:
                for m in Ml.get(g, []):
                    v = J.D({m: 1}, u)
                    if v:
                        vecs.append(v)
            self._DM[g] = vecs
        return self._DM[g]

    def QM(self, g):
        return [self.th.Q({m: 1}) for m in self.Mg(g)]


def rk3(vecs):
    if not vecs:
        return (0, 0, 0)
    return sla.ranks(vecs)


def cohomology(piece, gs, verbose=True):
    """returns dict g -> dict(dimM, rkD, rkQbar, H) for g in gs (needs g-1, g, g+1)."""
    out = {}
    rD = {}
    rQD = {}
    need = set()
    for g in gs:
        need.update([g - 1, g, g + 1])
    for g in sorted(need):
        rD[g] = rk3(piece.DM(g))
    for g in sorted(need):
        if g + 1 in need:
            rQD[g] = rk3([v for v in piece.QM(g) if v] + piece.DM(g + 1))
    for g in gs:
        dimM = len(piece.Mg(g))
        rkQ = tuple(rQD[g][i] - rD[g + 1][i] for i in range(3)) if g in rQD else None
        rkQm = tuple(rQD[g - 1][i] - rD[g][i] for i in range(3)) if (g - 1) in rQD else (0, 0, 0)
        H = tuple(dimM - rD[g][i] - rkQ[i] - rkQm[i] for i in range(3))
        out[g] = dict(dimM=dimM, rkD=rD[g], rkQbar=rkQ, rkQbar_prev=rkQm, H=H)
    return out
