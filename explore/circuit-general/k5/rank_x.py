"""Rank analysis of the K6 (completed K5) relation system at level q, modular arithmetic.
Unknown columns: every non-(T1 or T2) class; T1, T2 classes are boundary data.
Relations: (X) vertex IBP (all ordered (x, x1, c)); optionally (DS) links from level q-1 with level q-1 known.
usage: python3 rank_x.py Q [ds]"""
import sys
import time
import pickle
import os

import k6
import linsys

HERE = os.path.dirname(os.path.abspath(__file__))


def load_classes(q):
    path = os.path.join(HERE, 'data', 'classes_q%d.pkl' % q)
    if os.path.exists(path):
        return pickle.load(open(path, 'rb'))
    cl = k6.classes(q)
    os.makedirs(os.path.join(HERE, 'data'), exist_ok=True)
    pickle.dump(cl, open(path, 'wb'))
    return cl


def x_rows(q, cl, cache=None):
    rows = []
    for b in cl:
        for items in k6.x_relations(b, q):
            terms = [it[1] for it in items]
            cs = k6.canon_batch(terms, q)
            r = {}
            for (coef, _), c in zip(items, cs):
                r[c] = r.get(c, 0) + coef
            r = {k: v for k, v in r.items() if v}
            if r:
                rows.append(r)
    return rows


def ds_rows(q, cl_lower):
    """rows at level q from DS with level q-1 class b: J_{q-1}(b) - sum ... = 0; key ('L', b) for level q-1"""
    rows = []
    for b in cl_lower:
        for c in range(6):
            items = k6.ds_relation(b, q, c)
            if items is None:
                continue
            terms = [t for _, t in items]
            cs = k6.canon_batch(terms, q)
            r = {('L', b): -1}
            for (coef, _), cc in zip(items, cs):
                r[cc] = r.get(cc, 0) + coef
            rows.append({k: v for k, v in r.items() if v})
    return rows


if __name__ == '__main__':
    q = int(sys.argv[1])
    use_ds = 'ds' in sys.argv[2:]
    t = time.time()
    cl = load_classes(q)
    keys = sorted(cl)
    kinds = {b: k6.kind(b, q) for b in keys}
    dec = {b: bool(k6.k5_decompletions(b, q)) for b in keys}
    print('q=%d classes %d: T1 %d, T2 %d, G %d (K5-decompletable G %d)' % (
        q, len(keys), sum(1 for b in keys if kinds[b] == 'T1'), sum(1 for b in keys if kinds[b] == 'T2'),
        sum(1 for b in keys if kinds[b] == 'G'), sum(1 for b in keys if kinds[b] == 'G' and dec[b])), flush=True)
    rows = x_rows(q, keys)
    print('X rows %d (%.1fs)' % (len(rows), time.time() - t), flush=True)
    if use_ds:
        cl1 = load_classes(q - 1)
        rds = ds_rows(q, sorted(cl1))
        print('DS rows %d' % len(rds), flush=True)
        rows += rds
    # columns
    allk = set()
    for r in rows:
        allk |= set(r)
    missing = [k for k in allk if not (isinstance(k, tuple) and len(k) == 2 and k[0] == 'L') and k not in cl]
    assert not missing, missing[:3]
    G = [b for b in keys if kinds[b] == 'G']
    Gd = [b for b in G if dec[b]]
    Ga = [b for b in G if not dec[b]]
    T = [b for b in keys if kinds[b] != 'G']
    lower = sorted(k for k in allk if isinstance(k, tuple) and len(k) == 2 and k[0] == 'L')
    col_order = Gd + Ga + T + lower
    cidx = {k: i for i, k in enumerate(col_order)}
    mrows = [{cidx[k]: linsys.to_mod(v) for k, v in r.items()} for r in rows]
    order = list(range(len(col_order)))
    piv = linsys.rref_mod(mrows, len(col_order), col_order=order)
    print('rank %d of %d columns' % (len(piv), len(col_order)), flush=True)
    # (a) homogeneous: all classes unknown (T and G), lower known
    nunk_all = len(Gd) + len(Ga) + len(T)
    det, _ = linsys.determined(dict((c, dict(r)) for c, r in piv.items()), range(nunk_all), order)
    print('(a) everything at level q unknown: determined %d of %d; nullity (level-q) %d' % (
        len(det), nunk_all, nunk_all - sum(1 for c in piv if c < nunk_all)), flush=True)
    # (b) T1, T2 known
    det, _ = linsys.determined(dict((c, dict(r)) for c, r in piv.items()), range(len(Gd) + len(Ga)), order)
    nd = sum(1 for c in det if c < len(Gd))
    print('(b) T1,T2 known: generic K5-decompletable determined %d of %d; auxiliary determined %d of %d' % (
        nd, len(Gd), sum(1 for c in det if c >= len(Gd)), len(Ga)), flush=True)
    und = [col_order[c] for c in range(len(Gd)) if c not in det]
    for b in und[:20]:
        print('   undetermined:', b, 'decompletions at', k6.k5_decompletions(b, q))
    print('time %.1fs' % (time.time() - t))
