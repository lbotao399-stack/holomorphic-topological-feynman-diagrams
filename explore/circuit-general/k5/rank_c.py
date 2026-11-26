"""Rank analysis of the level-q K6 relation system (modular, C: modrank.c).
Modes:
  A  (X) only, every level-q class unknown: rank, nullity
  B  (X) only, T1 and T2 classes known: generic classes determined?
  B1 (X) only, T1 known only;  B2 (X) only, T2 known only
  Bd (X) only, only the K5-decompletable T1 and T2 classes known (auxiliary T1/T2 classes unknown)
  C  (X) + (DS from level q-1, level q-1 known), every level-q class unknown
usage: python3 rank_c.py Q [modes]"""
import os
import subprocess
import sys
import time

import k6
from solve_level import classes, class_file, SCR, HERE

P = 2147483647


def build():
    exe = os.path.join(SCR, 'modrank')
    src = os.path.join(HERE, 'modrank.c')
    if not os.path.exists(exe) or os.path.getmtime(exe) < os.path.getmtime(src):
        subprocess.check_call(['cc', '-O2', '-o', exe, src])
    return exe


def run(q, unknown_mask, nod, rf, nl):
    exe = build()
    n = len(unknown_mask)
    U = [j for j in range(n) if unknown_mask[j]]
    K = [j for j in range(n) if not unknown_mask[j]]
    perm = [0] * n
    for i, j in enumerate(U + K):
        perm[j] = i
    pf = os.path.join(SCR, 'perm_q%d.txt' % q)
    open(pf, 'w').write('\n'.join(map(str, perm)) + '\n')
    lf = os.path.join(SCR, 'low0_q%d.txt' % q)
    open(lf, 'w').write('%d 0\n' % nl)
    env = dict(os.environ, MS_PERM=pf, MS_NU=str(len(U)))
    if nod:
        env['MS_NOD'] = '1'
    out = subprocess.check_output([exe, str(P), str(n), rf, lf], env=env).decode().splitlines()
    head = out[0].split()
    rank = int(head[1])
    det = {}
    for line in out[1:]:
        c, d = map(int, line.split())
        det[U[c]] = d
    return rank, det


if __name__ == '__main__':
    q = int(sys.argv[1])
    modes = sys.argv[2].split(',') if len(sys.argv) > 2 else ['A', 'B', 'B1', 'B2', 'Bd', 'C']
    t0 = time.time()
    cls = classes(q)
    nl = len(classes(q - 1))
    rf = os.path.join(SCR, 'rows%d.txt' % q)
    if not os.path.exists(rf):
        from solve_level import rows_file
        rf = rows_file(q)
    n = len(cls)
    kind = [k6.kind(b, q) for b in cls]
    dec = [bool(k6.k5_decompletions(b, q)) for b in cls]
    nG = sum(1 for k in kind if k == 'G')
    nGd = sum(1 for j in range(n) if kind[j] == 'G' and dec[j])
    print('q=%d: classes %d (T1 %d, T2 %d, G %d of which K5-decompletable %d)' % (
        q, n, kind.count('T1'), kind.count('T2'), nG, nGd), flush=True)
    for m in modes:
        if m == 'A':
            mask, nod = [1] * n, True
        elif m == 'B':
            mask, nod = [1 if k == 'G' else 0 for k in kind], True
        elif m == 'B1':
            mask, nod = [0 if k == 'T1' else 1 for k in kind], True
        elif m == 'Bd':
            mask, nod = [0 if (kind[j] != 'G' and dec[j]) else 1 for j in range(n)], True
        elif m == 'B2':
            mask, nod = [0 if k == 'T2' else 1 for k in kind], True
        else:
            mask, nod = [1] * n, False
        rank, det = run(q, mask, nod, rf, nl)
        U = [j for j in range(n) if mask[j]]
        nd = sum(det.values())
        ndG = sum(1 for j in U if det[j] and kind[j] == 'G' and dec[j])
        nGU = sum(1 for j in U if kind[j] == 'G' and dec[j])
        print('  mode %-2s: unknown %d, rank %d, determined %d of %d; generic K5-decompletable determined %d of %d '
              '(%.0fs)' % (m, len(U), rank, nd, len(U), ndG, nGU, time.time() - t0), flush=True)
