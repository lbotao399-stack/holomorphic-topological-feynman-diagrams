"""Exact sparse linear algebra over Q (Fractions) with combination tracking."""
from fractions import Fraction as Fr


class Echelon:
    def __init__(self):
        self.rows = []      # (pivot, vec, combo)
        self.piv = {}
        self.ngen = 0

    def reduce(self, v, combo=None):
        v = dict(v)
        combo = dict(combo) if combo else {}
        for (p, r, c) in self.rows:
            if p in v:
                f = v[p] / r[p]
                for k, x in r.items():
                    y = v.get(k, 0) - f * x
                    if y == 0:
                        v.pop(k, None)
                    else:
                        v[k] = y
                for k, x in c.items():
                    y = combo.get(k, 0) - f * x
                    if y == 0:
                        combo.pop(k, None)
                    else:
                        combo[k] = y
        return v, combo

    def add(self, v, tag=None):
        """returns (True, None) if independent, else (False, relation combo) with
        v - sum combo[g] gen_g = 0 encoded as combo including tag:1."""
        if tag is None:
            tag = self.ngen
        self.ngen += 1
        res, combo = self.reduce(v, {tag: Fr(1)})
        if not res:
            return False, combo
        p = next(iter(res))
        self.rows.append((p, res, combo))
        return True, None

    def rank(self):
        return len(self.rows)

    def coords(self, v):
        """if v in span: dict gen->coef with v = sum coef gen; else None"""
        res, combo = self.reduce(v, {})
        if res:
            return None
        return {k: -x for k, x in combo.items()}


def rank_of(vecs):
    E = Echelon()
    for v in vecs:
        E.add(v)
    return E.rank()


def kernel(vecs):
    """relations sum c_i vecs_i = 0; returns list of dicts index->coef (basis of kernel)."""
    E = Echelon()
    rels = []
    for i, v in enumerate(vecs):
        ok, rel = E.add(v, i)
        if not ok:
            rels.append(rel)
    return rels
