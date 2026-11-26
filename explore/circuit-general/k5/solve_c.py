"""Level q of the K6 family from level q-1 by (X) + (DS), modular C solver (modsolve.c) over several primes,
CRT + rational reconstruction, verification with an independent prime and on all rows (exact, Fractions).
Undetermined classes get symbolic labels 'u<q>_<n>'.
With env K5FAST=1 (modsolve2.c): the first prime processes all rows and records the pivot rows; the other
primes process only those rows (the exact all-row check below is unchanged).
With env K5SUB=k: the C solver sees all DS rows and every k-th X row (the exact check uses all rows).
usage: python3 solve_c.py Q [nprimes]"""
import os
import subprocess
import sys
import time
from fractions import Fraction as F

import k6
from solve_level import classes, class_file, rows_file, read_rows, load_values, save_values, SCR, HERE

PRIMES = [2147483647, 2147483629, 2147483587, 2147483579, 2147483563, 2147483549, 2147483543, 2147483497,
          2147483489, 2147483477, 2147483423, 2147483399]


def build():
    name = 'modsolve2' if os.environ.get('K5FAST') else 'modsolve'
    exe = os.path.join(SCR, name)
    src = os.path.join(HERE, name + '.c')
    if not os.path.exists(exe) or os.path.getmtime(exe) < os.path.getmtime(src):
        subprocess.check_call(['cc', '-O2', '-o', exe, src])
    return exe


def tomod(x, p):
    return x.numerator % p * pow(x.denominator % p, p - 2, p) % p


def ratrec(a, m):
    """rational reconstruction of a mod m (|num|, den < sqrt(m/2))"""
    import math
    bound = math.isqrt(m // 2)
    r0, r1 = m, a % m
    s0, s1 = 0, 1
    while r1 > bound:
        qq = r0 // r1
        r0, r1 = r1, r0 - qq * r1
        s0, s1 = s1, s0 - qq * s1
    if s1 == 0 or abs(s1) > bound:
        return None
    return F(r1, s1) if s1 > 0 else F(-r1, -s1)


def run(q, nprimes=4, verbose=True):
    t0 = time.time()
    exe = build()
    cls = classes(q)
    low = classes(q - 1)
    lowvals = load_values(q - 1)
    labs = sorted({l for v in lowvals.values() for l in v}, key=lambda s: (s[0] == 'u', len(s), s))
    rf = rows_file(q)
    rfs = rf
    sub = int(os.environ.get('K5SUB', '1'))
    if sub > 1:
        rfs = rf + '.sub%d' % sub
        with open(rfs, 'w') as fo:
            for i, line in enumerate(open(rf)):
                if line[0] == 'D' or i % sub == 0:
                    fo.write(line)
    n = len(cls)
    res = []
    for p in PRIMES[:nprimes + 1]:
        lf = os.path.join(SCR, 'low_q%d_%d.txt' % (q, p))
        with open(lf, 'w') as f:
            f.write('%d %d\n' % (len(low), len(labs)))
            for b in low:
                if b in lowvals:
                    v = lowvals[b]
                    f.write(' '.join(str(tomod(F(v.get(l, 0)), p)) for l in labs) + '\n')
                else:
                    f.write(' '.join(['-1'] * len(labs)) + '\n')
        env = dict(os.environ)
        pr = os.path.join(SCR, 'pivrows_q%d.txt' % q)
        if os.environ.get('K5FAST'):
            if p == PRIMES[0]:
                env['MS_ROWS_OUT'] = pr
            else:
                env['MS_ROWS_IN'] = pr
        out = subprocess.check_output([exe, str(p), str(n), rfs, lf], env=env).decode().splitlines()
        head = out[0]
        det = {}
        for line in out[1:]:
            x = list(map(int, line.split()))
            if x[1]:
                det[x[0]] = x[2:]
        res.append((p, head, det))
        if verbose:
            print('q=%d p=%d: %s, determined %d of %d (%.1fs)' % (q, p, head, len(det), n, time.time() - t0), flush=True)
    dets = set(res[0][2])
    for p, h, d in res[1:]:
        assert set(d) == dets
    # CRT over the first nprimes, check with the last
    vals = {}
    nfail = 0
    for j in dets:
        v = {}
        for li, lab in enumerate(labs):
            a, m = 0, 1
            for p, h, d in res[:-1]:
                x = d[j][li]
                # CRT combine a mod m with x mod p
                t = (x - a) * pow(m, p - 2, p) % p
                a, m = a + m * t, m * p
            r = ratrec(a, m)
            if r is None:
                nfail += 1
                r = F(0)
            pl, hl, dl = res[-1]
            if tomod(r, pl) != dl[j][li]:
                nfail += 1
            if r:
                v[lab] = r
        vals[j] = v
    und = [j for j in range(n) if j not in dets]
    for k, j in enumerate(und):
        vals[j] = {'u%d_%d' % (q, k): F(1)}
    if verbose:
        print('q=%d: reconstruction failures %d; undetermined %d %s (%.1fs)' % (
            q, nfail, len(und), [cls[j] for j in und][:10], time.time() - t0), flush=True)
    # exact check on all rows (streamed)
    lowidx = {i: b for i, b in enumerate(low)}
    nbad = nchk = 0
    for line in open(rf):
        p_ = line.split()
        r = {}
        if p_[0] == 'X':
            items = p_[1:]
        else:
            r[('L', int(p_[2][1:]))] = F(int(p_[1]))
            items = p_[3:]
        for i in range(0, len(items), 2):
            j = int(items[i + 1])
            r[j] = r.get(j, 0) + int(items[i])
        acc = {}
        ok = True
        for k, c in r.items():
            if not c:
                continue
            if isinstance(k, tuple):
                b = lowidx[k[1]]
                if b not in lowvals:
                    ok = False
                    break
                vv = lowvals[b]
            else:
                vv = vals[k]
            for lab, x in vv.items():
                acc[lab] = acc.get(lab, 0) + c * x
        if not ok:
            continue
        nchk += 1
        if any(acc.values()):
            nbad += 1
    if verbose:
        print('q=%d: exact check on %d rows: %d fail (%.1fs)' % (q, nchk, nbad, time.time() - t0), flush=True)
    save_values(q, cls, vals, os.environ.get('K5TAG', 'tower'))
    return cls, vals, und, nbad, nfail


if __name__ == '__main__':
    q = int(sys.argv[1])
    npr = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    cls, vals, und, nbad, nfail = run(q, npr)
    labs = set()
    stat = {}
    for j, b in enumerate(cls):
        labs |= set(vals[j])
        kd = k6.kind(b, q)
        dec = bool(k6.k5_decompletions(b, q))
        key = kd + ('d' if dec else 'a')
        stat[key] = stat.get(key, 0) + 1
        if any(l.startswith('u') for l in vals[j]) and dec:
            print('   K5-decompletable class depends on an undetermined label:', b, vals[j])
    print('q=%d kinds %s; labels %s' % (q, stat, sorted(labs, key=lambda s: (s[0] == 'u', len(s), s))))
