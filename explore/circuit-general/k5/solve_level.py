"""Solve level q of the K6 family from the level-(q-1) values with (X) + (DS) rows (C generator xrows.c).
Undetermined level-q classes are carried as new symbolic labels 'u<q>_<n>' (J = u<q>_<n>).
Lower classes without a value are dropped from the DS rows that contain them.
Exact: modular selection of an independent row set (on the unknown columns), Fraction elimination of
the selected rows, check of every row on the solution.
usage: python3 solve_level.py Q"""
import json
import os
import subprocess
import sys
import time
from fractions import Fraction as F

import k6
import linsys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')
SCR = os.environ.get('K5SCR', '/private/tmp/claude-501/-Users-libotao-Desktop/07c118b5-2480-4b3c-899c-a49c36ac4d00/scratchpad')


def build_bins():
    for name in ('enum6', 'xrows'):
        exe = os.path.join(SCR, name)
        src = os.path.join(HERE, name + '.c')
        if not os.path.exists(exe) or os.path.getmtime(exe) < os.path.getmtime(src):
            subprocess.check_call(['cc', '-O2', '-o', exe, src])


def class_file(q):
    """sorted canonical class list (text) produced by enum6"""
    path = os.path.join(DATA, 'classes_q%d.txt' % q)
    if not os.path.exists(path):
        build_bins()
        out = subprocess.check_output([os.path.join(SCR, 'enum6'), str(q)]).decode()
        rows = sorted(tuple(int(x) for x in l.split()[:15]) for l in out.splitlines() if l.strip())
        os.makedirs(DATA, exist_ok=True)
        with open(path, 'w') as f:
            for r in rows:
                f.write(' '.join(map(str, r)) + '\n')
    return path


def classes(q):
    return [tuple(int(x) for x in l.split()) for l in open(class_file(q))]


def rows_file(q):
    path = os.path.join(SCR, 'rows%d.txt' % q)
    build_bins()
    with open(path, 'w') as f:
        subprocess.check_call([os.path.join(SCR, 'xrows'), str(q), class_file(q), class_file(q - 1)], stdout=f)
    return path


def read_rows(path):
    rows = []
    for line in open(path):
        p = line.split()
        r = {}
        if p[0] == 'X':
            items = p[1:]
        else:
            r[('L', int(p[2][1:]))] = F(int(p[1]))
            items = p[3:]
        for i in range(0, len(items), 2):
            c = int(items[i])
            j = int(items[i + 1])
            r[j] = r.get(j, 0) + c
        r = {k: v for k, v in r.items() if v}
        if r:
            rows.append(r)
    return rows


def save_values(q, cls, vals, tag='tower'):
    path = os.path.join(DATA, 'J_q%d_%s.json' % (q, tag))
    json.dump([[list(cls[i]), {k: str(x) for k, x in v.items()}] for i, v in sorted(vals.items())], open(path, 'w'))


def load_values(q, tag='tower'):
    path = os.path.join(DATA, 'J_q%d_%s.json' % (q, tag))
    if not os.path.exists(path):
        return None
    return {tuple(b): {k: F(x) for k, x in v.items()} for b, v in json.load(open(path))}


def select_rows(rows, unknown_set, p=linsys.P, order_key=None):
    piv = {}
    keep = []
    for ri, r in enumerate(rows):
        m = {}
        for k, v in r.items():
            if k in unknown_set:
                x = linsys.to_mod(v, p)
                if x:
                    m[k] = x
        while m:
            c = min(m, key=order_key) if order_key else min(m)
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


def solve(q, verbose=True):
    t0 = time.time()
    cls = classes(q)
    low = classes(q - 1)
    lowvals = load_values(q - 1)
    assert lowvals is not None, 'level %d values missing' % (q - 1)
    rows = read_rows(rows_file(q))
    nX = sum(1 for r in rows if not any(isinstance(k, tuple) for k in r))
    known = {}
    for i, b in enumerate(low):
        if b in lowvals:
            known[('L', i)] = lowvals[b]
    rows = [r for r in rows if all((not isinstance(k, tuple)) or k in known for k in r)]
    unknown = list(range(len(cls)))
    if verbose:
        print('q=%d: classes %d, rows %d (X %d), lower classes %d (with value %d)  (%.1fs)' % (
            q, len(cls), len(rows), nX, len(low), len(known), time.time() - t0), flush=True)
    keep, rk = select_rows(rows, set(unknown))
    if verbose:
        print('q=%d: modular rank %d of %d unknowns (%.1fs)' % (q, rk, len(unknown), time.time() - t0), flush=True)
    sel = [rows[i] for i in keep]
    sol, bad, rank = linsys.solve_exact(sel, unknown, known)
    und = [j for j in unknown if j not in sol]
    for n, j in enumerate(und):
        sol[j] = {'u%d_%d' % (q, n): F(1)}
    full = dict(known)
    full.update(sol)
    nbad = 0
    for r in rows:
        acc = {}
        for k, c in r.items():
            for lab, x in full[k].items():
                acc[lab] = acc.get(lab, 0) + c * x
        if any(acc.values()):
            nbad += 1
    if verbose:
        print('q=%d: exact rank %d, undetermined %d %s; all %d rows checked on the solution: %d fail (%.1fs)' % (
            q, rank, len(und), [cls[j] for j in und], len(rows), nbad, time.time() - t0), flush=True)
    vals = {j: v for j, v in sol.items()}
    save_values(q, cls, vals)
    return cls, vals, und, nbad


if __name__ == '__main__':
    q = int(sys.argv[1])
    cls, vals, und, nbad = solve(q)
    labs = set()
    for v in vals.values():
        labs |= set(v)
    stat = {}
    for j, b in enumerate(cls):
        kd = k6.kind(b, q)
        dec = bool(k6.k5_decompletions(b, q))
        key = kd + ('d' if dec else 'a')
        stat[key] = stat.get(key, 0) + 1
        if any(l.startswith('u') for l in vals[j]) and kd == 'G' and dec:
            print('   generic K5 class depends on an undetermined label:', b, vals[j])
    print('q=%d kinds %s; labels %s' % (q, stat, sorted(labs, key=lambda s: (s[0] == 'u', len(s), s))))
