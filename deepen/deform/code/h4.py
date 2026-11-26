"""H^4 of g[[z1,z2]] mod translations, invariant cochains, weight (2,2), formal gl_N and sl_N
(trace basis; faithful for N >= 5), exact ranks with certificates; classes of F_133, F_223.
Saves the H^4 data (bases, dual functionals) to h4_data.pkl."""
import sys, time, pickle
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/deform/code')
import ce, la, spaces

t0 = time.time()
obs = spaces.obstruction_cochains()
for pqr in obs:
    print("theta%s: (60/pi) K / ([p1p2][p3p4]) =" % (pqr,), obs[pqr]['Kcoef'])

DATA = {}
for sl in (False, True):
    name = "sl_N" if sl else "gl_N"
    print("\n==== formal %s, weight (2,2) ====" % name, flush=True)
    C3, keys3 = spaces.basis(3, (2, 2), sl)
    C4, keys4 = spaces.basis(4, (2, 2), sl)
    print("dim C^3_inv(2,2) = %d, dim C^4_inv(2,2) = %d  (ambient trace coordinates: %d, %d)" % (
        len(C3), len(C4), len(keys3), len(keys4)))
    # B^4
    dC3 = [spaces.const_part(ce.d_ce(spaces.to_np(c), 3, sl)) for c in C3]
    _, MB = la.matrix_from(dC3, keys4)
    cB = la.certified_rank(MB)
    print("rank d: C^3 -> C^4 = %d  (minor det %d, %d exact kernel vectors checked, mod-p ranks %s)" % (
        cB['rank'], cB['det'], len(cB['kernel']), cB.get('ranks_modp')))
    # Z^4
    dC4 = [spaces.const_part(ce.d_ce(spaces.to_np(c), 4, sl)) for c in C4]
    keys5 = spaces.keys_of(5, (2, 2), sl)
    _, MZ = la.matrix_from(dC4, keys5)
    cZ = la.certified_rank(MZ)
    print("rank d: C^4 -> C^5 = %d  (minor det %d, %d exact kernel vectors checked, mod-p ranks %s; C^5 coordinates %d)" % (
        cZ['rank'], cZ['det'], len(cZ['kernel']), cZ.get('ranks_modp'), len(keys5)))
    Z4 = []
    for v in cZ['kernel']:
        z = {}
        for j, cj in enumerate(v):
            if cj:
                for k, x in C4[j].items():
                    z[k] = z.get(k, 0) + cj * x
        z = {k: x for k, x in z.items() if x}
        Z4.append(z)
    # check Z4 closed exactly
    allclosed = all(ce.d_ce(spaces.to_np(z), 4, sl) == {} for z in Z4)
    B4 = [dC3[j] for j in cB['pivcols']]
    print("dim Z^4 = %d (all exactly closed: %s), dim B^4 = %d, dim H^4 = %d" % (
        len(Z4), allclosed, len(B4), len(Z4) - len(B4)))
    # obstruction cochains
    Fs = {}
    for pqr in obs:
        F = obs[pqr]['F']
        if sl:
            F = ce.restrict_sl(F)
        pw = sorted({p for d in F.values() for p in d})
        parts = {p: {k: d[p] for k, d in F.items() if d.get(p)} for p in pw}
        Fs[pqr] = parts
        print("F_%d%d%d: N-powers %s; d F = 0 (formal): %s" % (pqr + (pw, ce.d_ce(F, 4, sl) == {})))
    # H^4 basis: start from B4, add classes greedily: first the N^1 part of F_223, then Z4 vectors
    cand = []
    for pqr in [(2, 2, 3), (1, 3, 3)]:
        for p, c in sorted(Fs[pqr].items(), key=lambda x: -x[0]):
            cand.append(("F_%d%d%d[N^%d]" % (pqr + (p,)), c))
    cand += [("z%d" % i, z) for i, z in enumerate(Z4)]
    chosen = []
    cur = list(B4)
    for nm, c in cand:
        _, M = la.matrix_from(cur + [c], keys4)
        r = la.certified_rank(M)['rank']
        if r > len(cur):
            cur.append(c)
            chosen.append((nm, c))
    print("H^4 basis (mod B^4):", [nm for nm, c in chosen])
    # dual functionals Phi on the ambient coordinates: Phi(B4) = 0, Phi(h_j) = delta
    Hb = [c for nm, c in chosen]
    G = B4 + Hb
    _, MG = la.matrix_from(G, keys4)
    cG = la.certified_rank(MG)
    assert cG['rank'] == len(G)
    rows = cG['pivrows']
    sub = [[MG[i][j] for j in range(len(G))] for i in rows]
    # solve Phi_R sub = [0 | I]: Phi_R = [0 I] sub^{-1}; compute row by row: sub^T phi^T = e
    Phi = []
    for h in range(len(Hb)):
        rhs = [Fr(0)] * len(B4) + [Fr(int(t == h)) for t in range(len(Hb))]
        subT = [[sub[i][j] for i in range(len(rows))] for j in range(len(G))]
        x = la.solve_square(subT, rhs)
        Phi.append({keys4[rows[i]]: x[i] for i in range(len(rows)) if x[i] != 0})

    def coords(c):
        return [sum(Fr(c.get(k, 0)) * v for k, v in ph.items()) for ph in Phi]
    # check Phi on basis
    okPhi = all(coords(b) == [0] * len(Hb) for b in B4) and all(
        coords(h) == [Fr(int(i == j)) for j in range(len(Hb))] for i, h in enumerate(Hb))
    print("dual functionals: Phi(B^4) = 0 and Phi(h_j) = delta_ij exactly:", okPhi)
    for pqr in obs:
        for p, c in sorted(Fs[pqr].items()):
            inB = la.in_span(MB, la.column(c, keys4), cB)[0]
            print("  [F_%d%d%d[N^%d]] coordinates %s ; exact (in B^4): %s" % (pqr + (p, coords(c), inB)))
    # class of F_133 + F_223 and of F_223 as polynomials in N
    for nm, pqrs in (("F_223", [(2, 2, 3)]), ("F_133", [(1, 3, 3)]), ("F_133+F_223", [(1, 3, 3), (2, 2, 3)])):
        cl = {}
        for pqr in pqrs:
            for p, c in Fs[pqr].items():
                for i, x in enumerate(coords(c)):
                    if x:
                        cl.setdefault(i, {})
                        cl[i][p] = cl[i].get(p, 0) + x
        print("  [%s](N) = %s" % (nm, {Hn: {p: str(v) for p, v in d.items() if v} for Hn, d in
                                       ((chosen[i][0], cl[i]) for i in cl)}))
    DATA[name] = dict(keys4=keys4, B4=B4, Z4=Z4, Hbasis=chosen, Phi=Phi, Fs=Fs, C4=C4, C3=C3)

pickle.dump(DATA, open((__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/deform/code/h4_data.pkl', 'wb'))
print("\ntotal %.0fs" % (time.time() - t0))
