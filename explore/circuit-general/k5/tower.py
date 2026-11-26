"""Level-by-level solution of the K6 (completed K5) family by (X) vertex IBP + (DS) dimension links.

Level 2: the single class (octahedron, = W4 at q=2), J = 20 z5 (engine and graphical functions).
Level 3: 15 classes, all T1/T2, values by bval (engine / graphical functions / per5).
Level q >= 4: unknown = all classes; rows = (X) at level q + (DS) from every level-(q-1) class and every c,
with level-(q-1) values known.  Optional boundary rows J = value for given classes.
Exact arithmetic: modular pre-selection of an independent row set, then Fraction elimination.
usage: python3 tower.py QMAX
"""
import json
import os
import pickle
import sys
import time
from fractions import Fraction as F

import k6
import linsys
from rank_x import load_classes, x_rows, ds_rows

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')


def vstr(v):
    order = sorted(v, key=lambda s: (s != '1', len(s), s))
    return ' + '.join('(%s)%s' % (v[k], '' if k == '1' else '*' + k) for k in order) or '0'


def save_values(q, vals, tag):
    path = os.path.join(DATA, 'J_q%d_%s.json' % (q, tag))
    json.dump([[list(b), {k: str(x) for k, x in v.items()}] for b, v in sorted(vals.items())], open(path, 'w'))


def load_values(q, tag):
    path = os.path.join(DATA, 'J_q%d_%s.json' % (q, tag))
    if not os.path.exists(path):
        return None
    return {tuple(b): {k: F(x) for k, x in v.items()} for b, v in json.load(open(path))}


def select_rows(rows, cols_unknown, p=linsys.P):
    """indices of a maximal set of rows independent on the unknown columns (modular), plus rank"""
    cidx = {k: i for i, k in enumerate(cols_unknown)}
    order = list(range(len(cols_unknown)))
    piv = {}
    keep = []
    for ri, r in enumerate(rows):
        m = {cidx[k]: linsys.to_mod(v, p) for k, v in r.items() if k in cidx}
        m = {c: v for c, v in m.items() if v}
        while m:
            c = min(m)
            if c in piv:
                f = m[c]
                for cc, vv in piv[c].items():
                    nv = (m.get(cc, 0) - f * vv) % p
                    if nv:
                        m[cc] = nv
                    else:
                        m.pop(cc, None)
            else:
                f = linsys.inv(m[c], p)
                piv[c] = {cc: vv * f % p for cc, vv in m.items()}
                keep.append(ri)
                break
    return keep, len(piv)


def solve_level(q, lower_vals, boundary=None, verbose=True):
    t0 = time.time()
    cl = load_classes(q)
    keys = sorted(cl)
    rows = x_rows(q, keys)
    nx = len(rows)
    if lower_vals is not None:
        rows += ds_rows(q, sorted(lower_vals))
    known = {}
    if lower_vals is not None:
        for b, v in lower_vals.items():
            known[('L', b)] = v
    if boundary:
        for b, v in boundary.items():
            known[b] = v
    unknown = [b for b in keys if b not in known]
    keep, rk = select_rows(rows, unknown)
    if verbose:
        print('q=%d: classes %d, X rows %d, DS rows %d, unknown %d, modular rank %d (%.1fs)' % (
            q, len(keys), nx, len(rows) - nx, len(unknown), rk, time.time() - t0), flush=True)
    sel = [rows[i] for i in keep]
    sol, bad, rank = linsys.solve_exact(sel, unknown, known)
    # consistency: all rows (not only selected) evaluated on the solution
    full = dict(known)
    full.update(sol)
    nbad = 0
    nchk = 0
    for r in rows:
        if all(k in full for k in r):
            nchk += 1
            acc = {}
            for k, c in r.items():
                for lab, x in full[k].items():
                    acc[lab] = acc.get(lab, 0) + c * x
            if any(acc.values()):
                nbad += 1
    if verbose:
        print('q=%d: exact rank %d, determined %d of %d, inconsistent selected rows %d; all-row check %d rows, %d fail '
              '(%.1fs)' % (q, rank, len(sol), len(unknown), len(bad), nchk, nbad, time.time() - t0), flush=True)
    return sol, len(unknown), nbad


if __name__ == '__main__':
    import bval
    qmax = int(sys.argv[1])
    os.makedirs(DATA, exist_ok=True)
    # level 2, 3: direct values
    for q in (2, 3):
        vals = load_values(q, 'direct')
        if vals is None:
            vals = {}
            for b in sorted(load_classes(q)):
                r = bval.value(b, q, methods=('E', 'S', 'P'))
                assert r is not None, b
                vals[b] = r[0]
            save_values(q, vals, 'direct')
        print('level %d: %d classes, all direct' % (q, len(vals)), flush=True)
    # check level-3 rows on direct values
    v2 = load_values(2, 'direct')
    v3 = load_values(3, 'direct')
    sol, nu, nbad = solve_level(3, v2, boundary=v3)
    print('level 3: X + DS(2->3) rows on direct values: %d failures' % nbad, flush=True)
    sol3, nu, nbad = solve_level(3, v2)
    agree = sum(1 for b in sol3 if sol3[b] == v3[b])
    print('level 3 from level 2 alone: determined %d of %d, agree with direct %d' % (len(sol3), nu, agree), flush=True)
    prev = v3
    for q in range(4, qmax + 1):
        vals = load_values(q, 'tower')
        if vals is None:
            sol, nu, nbad = solve_level(q, prev)
            if len(sol) < nu:
                print('level %d NOT closed: %d of %d' % (q, len(sol), nu), flush=True)
            vals = sol
            save_values(q, vals, 'tower')
        labs = set()
        for v in vals.values():
            labs |= set(v)
        print('level %d: %d values; labels %s' % (q, len(vals), sorted(labs, key=lambda s: (len(s), s))), flush=True)
        prev = vals
