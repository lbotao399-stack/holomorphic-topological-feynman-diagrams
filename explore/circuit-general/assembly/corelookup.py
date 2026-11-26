"""Vectorized exact core periods  I_q(G; alpha) = int prod t^(alpha-1) / U^q  for many alpha.

G: series skeleton (multigraph).  Structural series-parallel reduction (powers enter linearly:
series s_a + s_b, parallel s_a + s_b - q) down to a final core; normalization of cores/msys.py:
  P(G; a) = prod_e Gamma(a_e) * prod_P G_2q(s_a, s_b) / prod_e' Gamma(a'_e') * P(G_final; a'),
  G_2q(a, b) = Gamma(q-a) Gamma(q-b) Gamma(a+b-q) / (Gamma(a) Gamma(b) Gamma(2q-a-b)),
  parallel step convergent iff 0 < s_a, s_b < q, s_a + s_b > q.
Final core: banana (P = prod Gamma(q - s_i) / Gamma(q)), prism / K33 (orbit tables of periods.py,
expanded by automorphisms, sorted-key lookup), wheel W_n (cores/wheel.py, memoized on the
dihedral-canonical (spoke, rim)).
prepare(A) is prime independent; evaluate(prep, p) returns residues per label and the
divergence mask.
"""
import os
import sys
from fractions import Fraction as F
from itertools import permutations

import numpy as np

import periods
import gutil
import msys
import vac
import wheel

_WHEEL_MEMO = {}
_EXPANDED = {}


def _encode(A, base):
    key = np.zeros(A.shape[0], dtype=np.int64)
    for j in range(A.shape[1]):
        key = key * base + A[:, j].astype(np.int64)
    return key


def expanded_table(name, q):
    """(sorted full keys, orbit index per key, list of value dicts, base)"""
    k = (name, q)
    if k in _EXPANDED:
        return _EXPANDED[k]
    tab = periods.table(name, q)
    fe = periods.FAMILIES[name]
    autos = gutil.edge_automorphisms(tuple(fe))
    reps = sorted(tab)
    R = np.array(reps, dtype=np.int64)
    base = 4 * q + 2
    keys, idx = [], []
    for pmt in autos:
        keys.append(_encode(R[:, list(pmt)], base))
        idx.append(np.arange(len(reps), dtype=np.int64))
    keys = np.concatenate(keys)
    idx = np.concatenate(idx)
    order = np.argsort(keys, kind='stable')
    keys, idx = keys[order], idx[order]
    uk, first = np.unique(keys, return_index=True)
    out = (uk, idx[first], [tab[r] for r in reps], base)
    _EXPANDED[k] = out
    return out


def wheel_value(spoke, rim, q):
    n = len(spoke)
    cands = []
    for s in range(n):
        for rev in (False, True):
            if not rev:
                sp = tuple(spoke[(s + k) % n] for k in range(n))
                rm = tuple(rim[(s + k) % n] for k in range(n))
            else:
                sp = tuple(spoke[(s - k) % n] for k in range(n))
                rm = tuple(rim[(s - k - 1) % n] for k in range(n))
            cands.append((sp, rm))
    key = (min(cands), q)
    v = _WHEEL_MEMO.get(key)
    if v is None:
        sp, rm = key[0]
        val = wheel.integrate(list(sp), list(rm), q)
        v = {('1' if k == 0 else 'z%d' % k): x for k, x in val.items() if x}
        _WHEEL_MEMO[key] = v
    return v


class CoreLookup:
    def __init__(self, edges, q):
        self.edges0 = [tuple(e) for e in edges]
        self.q = q
        E = len(edges)
        # current edges: endpoints, power = mask . alpha - m q
        cur = [dict(ends=tuple(e), mask=np.eye(E, dtype=np.int64)[i], m=0) for i, e in enumerate(edges)]
        self.pops = []          # (mask_a, m_a, mask_b, m_b)
        while True:
            verts = sorted({x for c in cur for x in c['ends']})
            if len(verts) == 2:
                break
            deg = {v: 0 for v in verts}
            for c in cur:
                for x in c['ends']:
                    deg[x] += 1
            done = False
            for v in verts:
                if deg[v] == 2:
                    i, j = [k for k, c in enumerate(cur) if v in c['ends']]
                    a = cur[i]['ends'][0] if cur[i]['ends'][1] == v else cur[i]['ends'][1]
                    b = cur[j]['ends'][0] if cur[j]['ends'][1] == v else cur[j]['ends'][1]
                    if a == b:
                        raise ValueError('tadpole in skeleton reduction')
                    new = dict(ends=(a, b), mask=cur[i]['mask'] + cur[j]['mask'], m=cur[i]['m'] + cur[j]['m'])
                    cur = [c for k, c in enumerate(cur) if k not in (i, j)] + [new]
                    done = True
                    break
            if done:
                continue
            seen = {}
            for i, c in enumerate(cur):
                pr = tuple(sorted(c['ends']))
                if pr in seen:
                    j = seen[pr]
                    self.pops.append((cur[j]['mask'], cur[j]['m'], c['mask'], c['m']))
                    new = dict(ends=pr, mask=cur[j]['mask'] + c['mask'], m=cur[j]['m'] + c['m'] + 1)
                    cur = [d for k, d in enumerate(cur) if k not in (i, j)] + [new]
                    done = True
                    break
                seen[pr] = i
            if not done:
                break
        self.final = cur
        fedges = [c['ends'] for c in cur]
        self.fedges = fedges
        verts = sorted({x for e in fedges for x in e})
        if len(verts) == 2:
            self.kind = 'banana'
        else:
            deg = {v: 0 for v in verts}
            for u, v in fedges:
                deg[u] += 1
                deg[v] += 1
            if vac.wheel_form(fedges, list(range(len(fedges))), verts, deg) is not None:
                self.kind = 'wheel'
                self.verts, self.deg = verts, deg
            else:
                sk = msys.struct_key(fedges)
                self.kind = None
                for name, fe in periods.FAMILIES.items():
                    if msys.struct_key(fe) == sk:
                        self.kind = name
                        S = msys.MultiSystem([fe], q)
                        # family vector f[j] = final power at edge fam_take[j]
                        self.fam_take = S.map_to_family(fedges, list(range(len(fedges))), 0)
                        break
                if self.kind is None:
                    raise ValueError('unsupported final core %s' % (fedges,))
        self.Fmask = np.array([c['mask'] for c in cur], dtype=np.int64)    # Ef x E
        self.Fm = np.array([c['m'] for c in cur], dtype=np.int64)
        self.n_loops = gutil.cycle_rank(self.edges0)
        self.sub_masks = None

    def _final_powers(self, A):
        return A.astype(np.int64) @ self.Fmask.T - self.q * self.Fm[None, :]

    def prepare(self, A):
        """A: (M x E) int powers (all >= 1).  Prime-independent data."""
        q = self.q
        A = np.asarray(A, dtype=np.int64)
        M = A.shape[0]
        ok = np.ones(M, dtype=bool)
        for ma, mqa, mb, mqb in self.pops:
            sa = A @ ma - q * mqa
            sb = A @ mb - q * mqb
            ok &= (sa > 0) & (sa < q) & (sb > 0) & (sb < q) & (sa + sb > q)
        Af = self._final_powers(A)
        prep = dict(M=M, ok=ok)
        if self.kind == 'banana':
            k = Af.shape[1]
            ok &= ((Af > 0) & (Af < q)).all(axis=1) & (Af.sum(axis=1) == (k - 1) * q)
            prep['idx'] = None
        elif self.kind in ('prism', 'K33'):
            uk, oidx, vals, base = expanded_table(self.kind, q)
            Fv = Af[:, self.fam_take]
            good = ((Fv >= 1) & (Fv < base)).all(axis=1)
            key = _encode(np.where(good[:, None], Fv, 1), base)
            pos = np.searchsorted(uk, key)
            pos = np.minimum(pos, len(uk) - 1)
            found = good & (uk[pos] == key)
            ok &= found
            prep['idx'] = np.where(found, oidx[pos], -1)
            prep['vals'] = vals
        else:   # wheel: persistent key -> value index cache (sorted arrays + dict for new keys)
            if not hasattr(self, '_wkeys'):
                self._wkeys, self._wvals = {}, []
                self._wbase = 4 * q + 2
                self._wk_sorted = np.zeros(0, dtype=np.int64)
                self._wk_loc = np.zeros(0, dtype=np.int64)
            okw = ok & (Af >= 1).all(axis=1) & (Af < self._wbase).all(axis=1)
            keys = _encode(np.where(okw[:, None], Af, 1), self._wbase)
            uk, first, inv = np.unique(keys, return_index=True, return_inverse=True)
            loc = np.full(len(uk), -2, dtype=np.int64)
            if len(self._wk_sorted):
                pos = np.minimum(np.searchsorted(self._wk_sorted, uk), len(self._wk_sorted) - 1)
                hit = self._wk_sorted[pos] == uk
                loc[hit] = self._wk_loc[pos[hit]]
            new = np.nonzero(loc == -2)[0]
            for u in new:
                kk = int(uk[u])
                pw = [int(x) for x in Af[first[u]]]
                if not gutil.convergent(tuple(self.fedges), pw, q):
                    r = -1
                else:
                    sp, rm = vac.wheel_form(self.fedges, pw, self.verts, self.deg)
                    self._wvals.append(wheel_value(sp, rm, q))
                    r = len(self._wvals) - 1
                self._wkeys[kk] = r
                loc[u] = r
            if len(new):
                allk = np.concatenate([self._wk_sorted, uk[new]])
                alll = np.concatenate([self._wk_loc, loc[new]])
                o = np.argsort(allk, kind='stable')
                self._wk_sorted, self._wk_loc = allk[o], alll[o]
            vidx = loc[inv.ravel()]
            vidx = np.where(okw, vidx, -1)
            ok &= vidx >= 0
            prep['idx'] = vidx
            prep['vals'] = self._wvals
        prep['ok'] = ok
        if prep.get('idx') is not None:
            prep['idx'] = prep['idx'].astype(np.int32)
        return prep

    def gamma_args(self, A):
        q = self.q
        A = np.asarray(A, dtype=np.int64)
        gnum, gden = [A], []
        for ma, mqa, mb, mqb in self.pops:
            sa = A @ ma - q * mqa
            sb = A @ mb - q * mqb
            gnum += [(q - sa)[:, None], (q - sb)[:, None], (sa + sb - q)[:, None]]
            gden += [sa[:, None], sb[:, None], (2 * q - sa - sb)[:, None]]
        Af = self._final_powers(A)
        gden.append(Af)
        if self.kind == 'banana':
            gnum.append(q - Af)
        return np.concatenate(gnum, axis=1), np.concatenate(gden, axis=1)

    def evaluate(self, A, prep, p, fact=None):
        """dict label -> int64 residues (0 on divergent rows), and the ok mask; A as in prepare"""
        q = self.q
        M, ok = prep['M'], prep['ok']
        gnum, gden = self.gamma_args(A)
        top = int(max(gnum.max(initial=1), gden.max(initial=1), q)) + 2
        if fact is None or len(fact) < top:
            fact = [1]
            for i in range(1, top):
                fact.append(fact[-1] * i % p)
        fa = np.array(fact, dtype=np.int64)
        ifa = np.array([pow(int(x), p - 2, p) for x in fa], dtype=np.int64)
        g = np.ones(M, dtype=np.int64)
        gn = np.where(ok[:, None], gnum, 1)
        gd = np.where(ok[:, None], gden, 1)
        del gnum, gden
        for j in range(gn.shape[1]):
            g = g * fa[gn[:, j] - 1] % p
        for j in range(gd.shape[1]):
            g = g * ifa[gd[:, j] - 1] % p
        out = {}
        if self.kind == 'banana':
            g = g * ifa[q - 1] % p
            out['1'] = np.where(ok, g, 0)
            return out, ok
        vals = prep['vals']
        idx = prep['idx']
        labels = sorted({k for v in vals for k in v})
        sel = np.where(ok, idx, len(vals))
        for lab in labels:
            vm = value_residues(vals, lab, p)
            out[lab] = vm[sel] * g % p
        return out, ok


_VRES = {}


def value_residues(vals, lab, p):
    """residues of one label of a list of value dicts, with a trailing 0 (cached per list)"""
    k = (id(vals), lab, p)
    r = _VRES.get(k)
    if r is None or r[0] is not vals or len(r[1]) != len(vals) + 1:
        arr = np.array([periods_fmod(v.get(lab, 0), p) for v in vals] + [0], dtype=np.int64)
        if len(_VRES) > 200:
            _VRES.clear()
        r = (vals, arr)
        _VRES[k] = r
    return r[1]


_FMOD = {}


def periods_fmod(x, p):
    if not x:
        return 0
    k = (x, p)
    r = _FMOD.get(k)
    if r is None:
        x = F(x)
        r = x.numerator % p * pow(x.denominator % p, p - 2, p) % p
        if len(_FMOD) > 5000000:
            _FMOD.clear()
        _FMOD[k] = r
    return r
