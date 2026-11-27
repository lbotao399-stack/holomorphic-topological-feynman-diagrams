"""Order-n nilpotency sum_{a+b=n} Q_a Q_b O = 0 with explicit operators (H <= 2, or small cases).

Q_n = c^n sum_Gamma kappa_Gamma hatQ_Gamma; c drops out (all terms of order n carry c^n).  At each order the
unknown kappa_n (one per graph) solve the linear system  sum_Gamma kappa_Gamma (E_{0,Gamma} + E_{Gamma,0})
= - sum_{a,b>=1} kappa kappa' E_{a,b}  over all coefficients of the random test operators; kappa_1 := sigma/|Aut| with
sigma = +1 for the first graph (normalization).  usage:
  python3 nilp_explicit.py THEORY GAUGE MATTER NMAX NOPS SEED [JETMAX] [NLET]
"""
import sys
import random
import time
from fractions import Fraction as Fr

import polyops as po
import theories
from wick import Algebra
from explicit import graph_ops


def rank_solve(rows, rhs, P):
    """rows: list of lists (mod P); returns (rank, augmented rank, solution or None)."""
    M = [list(r) + [b] for r, b in zip(rows, rhs)]
    nc = len(rows[0]) if rows else 0
    rk = 0
    pivcols = []
    for c in range(nc):
        piv = next((i for i in range(rk, len(M)) if M[i][c] % P), None)
        if piv is None:
            continue
        M[rk], M[piv] = M[piv], M[rk]
        inv = pow(M[rk][c], P - 2, P)
        M[rk] = [x * inv % P for x in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][c] % P:
                fac = M[i][c]
                M[i] = [(x - fac * y) % P for x, y in zip(M[i], M[rk])]
        pivcols.append(c)
        rk += 1
    aug = rk + (1 if any(M[i][nc] % P for i in range(rk, len(M))) else 0)
    sol = None
    if rk == nc and aug == rk:
        sol = [M[i][nc] for i in range(nc)]
    return rk, aug, sol


def rand_ops(alg, rng, k, nlet, jetmax):
    out = []
    while len(out) < k:
        lets = [(rng.randrange(alg.N), tuple(rng.randint(0, jetmax) for _ in range(alg.H))) for _ in range(nlet)]
        s, m = po.canon(lets, alg.pi)
        if s:
            out.append({m: 1})
    return out


def main(argv):
    th, gname, matter = argv[0], argv[1], (None if argv[2] == '-' else argv[2])
    nmax, nops, seed = int(argv[3]), int(argv[4]), int(argv[5])
    argv = [a for a in argv if a != '--unicyclic'] + (['--unicyclic'] if '--unicyclic' in argv else [])
    jetmax = int(argv[6]) if len(argv) > 6 and argv[6] != '--unicyclic' else 1
    nlet = int(argv[7]) if len(argv) > 7 and argv[7] != '--unicyclic' else None
    H, sl = theories.build(th, gname, matter)
    alg = Algebra(sl, H)
    rng = random.Random(seed)
    P = po.P
    print('theory %s gauge %s matter %s H=%d N=%d nmax=%d ops=%d seed=%d jetmax=%d P=%d'
          % (th, gname, matter, H, alg.N, nmax, nops, seed, jetmax, P), flush=True)
    uni = '--unicyclic' in argv
    ops = {m: graph_ops(alg, m, unicyclic=uni) for m in range(1, nmax + 1)}
    for m in ops:
        print('  order %d graphs: %s' % (m, [(g.name, g.aut, g.edges) for g in ops[m]]), flush=True)
    Q0 = alg.q0
    kap = {}
    T0 = time.time()
    for n in range(1, nmax + 1):
        # test operators with n+1 .. n+2 letters
        tests = []
        for L in ([nlet] if nlet else [n + 1, n + 2]):
            got, tries = 0, 0
            while got < nops and tries < 400 * nops:
                tries += 1
                O = rand_ops(alg, rng, 1, L, jetmax)[0]
                if any(g.apply(O) for g in ops[n]):
                    tests.append(O)
                    got += 1
        rows, rhs = [], []
        names = [g.name for g in ops[n]]
        for O in tests:
            img0 = Q0(O)
            cols = []
            for g in ops[n]:
                e = po.psum(Q0(g.apply(O)), g.apply(img0))
                cols.append(e)
            r = {}
            for a in range(1, n):
                b = n - a
                for ga in ops[a]:
                    for gb in ops[b]:
                        e = ga.apply(gb.apply(O))
                        k = kap[ga.name] * kap[gb.name] % P
                        for mo, v in e.items():
                            po.padd(r, mo, -k * v)
            monos = set(r)
            for c in cols:
                monos |= set(c)
            for mo in monos:
                rows.append([c.get(mo, 0) for c in cols])
                rhs.append(r.get(mo, 0))
        if n == 1:
            # homogeneous: check {Q0, Q1} = 0 for each graph separately
            nz = sum(1 for row in rows if any(x % P for x in row))
            nontriv = sum(1 for O in tests if any(g.apply(O) for g in ops[1]))
            print('order 1: {Q0, hatQ_G} O = 0 on %d test operators: %s  (nonzero rows %d); hatQ_1 O != 0 on %d'
                  % (len(tests), nz == 0, nz, nontriv), flush=True)
            for g in ops[1]:
                kap[g.name] = po.modq(Fr(1, g.aut))
            continue
        rk, aug, sol = rank_solve(rows, rhs, P)
        print('order %d: %d equations, rank %d, augmented rank %d, unknowns %d' % (n, len(rows), rk, aug, len(names)),
              flush=True)
        if sol is None:
            print('  NO SOLUTION / underdetermined', flush=True)
            return
        for nm, s, g in zip(names, sol, ops[n]):
            kap[nm] = s
            print('  kappa[%s] = %s   |Aut| = %d   kappa*|Aut| = %s' % (nm, po.ratrec(s), g.aut,
                                                                       po.ratrec(s * g.aut % P)), flush=True)
        print('  time %.1fs' % (time.time() - T0), flush=True)


if __name__ == '__main__':
    main(sys.argv[1:])
