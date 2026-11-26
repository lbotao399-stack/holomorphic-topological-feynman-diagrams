"""Per-algebra computation: gl1..gl5, sl2..sl5.
(a) H^3 at weights (1,1), (2,0), (0,2);  (b) classes [(1/2){J,J}] in H^4(2,2) in a basis starting with
[F_223] (if nonzero), membership of [F_223] in the image of the quadratic map (Groebner basis over Q,
Nullstellensatz with the Rabinowitsch variable);  (c) dim H^4(2,2).
Evaluation at N: formal trace coordinates when faithful (N >= k), else explicit.Evaluator."""
import sys, time, pickle
from fractions import Fraction as Fr
import sympy as sp
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/deform/code')
import ce, la, spaces, explicit

t00 = time.time()
obs = spaces.obstruction_cochains()
F223_formal = obs[(2, 2, 3)]['F']
F133_formal = obs[(1, 3, 3)]['F']

ALGS = [(1, False), (2, False), (3, False), (4, False), (5, False), (2, True), (3, True), (4, True), (5, True)]
if len(sys.argv) > 1:
    ALGS = [a for a in ALGS if ("sl%d" if a[1] else "gl%d") % a[0] in sys.argv[1:]]
RESULTS = {}


class Coord:
    """faithful coordinates of k-cochains at (N, sl)."""
    def __init__(self, N, sl, k):
        self.N, self.sl, self.k = N, sl, k
        self.ev = None
        if N < k:
            if k == 5 and N == 4:
                ev = explicit.Evaluator(N, sl, k, mode='samples', nsamp=(160 if not sl else 70), seed=11)
                rel = explicit.antisym_relation(k, sl)
                okrel = explicit.check_relation_on_samples(rel, N, sl, k)
                ntau = len(ev.taus)
                ev.exact = okrel and ev.rank == ntau - 1
                print("    antisymmetrizer identity (k=5) holds at N=4 on 6 fresh samples: %s; rank %d = %d - 1: exact %s" % (
                    okrel, ev.rank, ntau, ev.exact))
                assert ev.exact
            else:
                ev = explicit.Evaluator(N, sl, k, mode='tensor')
            self.ev = ev

    def vec(self, c):
        if self.ev is not None:
            return self.ev.vec(c)
        out = {}
        for (t, m), v in c.items():
            if self.sl and ce.has_singleton(t):
                continue
            x = ce.np_eval(v, self.N) if isinstance(v, dict) else Fr(v)
            if x != 0:
                out[(t, m)] = x
        return out


def basis_from(cochains, coord):
    vecs = [coord.vec(c) for c in cochains]
    keys, M = la.matrix_from(vecs)
    cert = la.certified_rank(M) if vecs and keys else dict(rank=0, pivcols=[], kernel=[])
    return [cochains[j] for j in cert['pivcols']], cert


def combo(cs, coefs):
    out = {}
    for c, x in zip(cs, coefs):
        if x:
            for k, d in c.items():
                ce.add_term(out, k, d if isinstance(d, dict) else {0: Fr(d)}, x)
    return out


def npl(c):
    return {k: (d if isinstance(d, dict) else {0: Fr(d)}) for k, d in c.items()}


for (N, sl) in ALGS:
    name = ("sl%d" if sl else "gl%d") % N
    t0 = time.time()
    print("\n================ %s ================" % name, flush=True)
    co3, co4, co5 = Coord(N, sl, 3), Coord(N, sl, 4), Coord(N, sl, 5)
    R = {}
    # ---------------- (a) H^3
    Z3 = {}
    for deg in [(1, 1), (2, 0), (0, 2)]:
        sp3 = [npl(c) for c in spaces.spanning(3, deg, sl)]
        C3, _ = basis_from(sp3, co3) if sp3 else ([], None)
        sp2 = spaces.spanning(2, deg, sl)
        C2 = [c for c in sp2 if co3 and co3.vec(c)] if False else sp2
        dC3 = [ce.d_ce(c, 3, sl) for c in C3]
        vecs = [co4.vec(c) for c in dC3]
        keys, M = la.matrix_from(vecs)
        if C3 and keys:
            cert = la.certified_rank(M)
            ker = cert['kernel']
            rk = cert['rank']
        else:
            ker = [[Fr(int(i == j)) for i in range(len(C3))] for j in range(len(C3))]
            rk = 0
        Z = [combo(C3, v) for v in ker]
        Z3[deg] = Z
        print("  (a) weight %s: dim C^2 = %d, dim C^3 = %d, rank d3 = %d, dim H^3 = dim Z^3 = %d" % (
            deg, len(C2), len(C3), rk, len(Z)), flush=True)
    R['dimH3'] = {deg: len(Z3[deg]) for deg in Z3}
    # GL(2) classes D, S1, S0 and their status at this N
    u = ce.bracket_poly(0, 1, 3)
    Ts = {"D": {((0, 1, 2),): {0: Fr(1)}, ((0, 2, 1),): {0: Fr(1)}},
          "S1": {((0,), (1, 2)): {0: Fr(1)}, ((0, 2), (1,)): {0: Fr(1)}, ((0, 1), (2,)): {0: Fr(1)}},
          "S0": {((0,), (1,), (2,)): {0: Fr(1)}}}
    if sl:
        Ts = {"D": Ts["D"]}
    GLJ = {nm: ce.canon_raw(ce.raw_from(u, T), 3, sl=sl) for nm, T in Ts.items()}
    # express a basis of Z^3(1,1) and check which GL(2) J's are nonzero / independent at this N
    Zb = Z3[(1, 1)]
    allJ = list(GLJ.values())
    keysZ, MZ = la.matrix_from([co3.vec(c) for c in Zb + allJ]) if (Zb or allJ) else ([], [])
    if keysZ:
        cz = la.certified_rank([row[:len(Zb)] for row in MZ]) if Zb else dict(rank=0)
        cJ = la.certified_rank([row[len(Zb):] for row in MZ])
        call = la.certified_rank(MZ)
        print("  GL(2) J's %s: rank %d; rank Z^3(1,1) %d; rank of both %d (J's span Z^3(1,1): %s)" % (
            list(GLJ), cJ['rank'], cz['rank'], call['rank'], call['rank'] == cJ['rank'] == cz['rank']), flush=True)
        R['GLJ_span'] = (cJ['rank'], cz['rank'], call['rank'])
    else:
        print("  Z^3(1,1) = 0 and the GL(2) J's vanish at this N")
        R['GLJ_span'] = (0, 0, 0)
    # use as basis of H^3(1,1) the independent GL(2) J's when they span, else Zb
    if keysZ and call['rank'] == cJ['rank'] == cz['rank']:
        Jnames = [list(GLJ)[j] for j in cJ['pivcols']]
        Jbasis = [GLJ[n_] for n_ in Jnames]
    else:
        Jnames = ["z%d" % i for i in range(len(Zb))]
        Jbasis = Zb
    print("  basis of H^3(1,1):", Jnames)
    # ---------------- (c) H^4(2,2)
    sp4 = [npl(c) for c in spaces.spanning(4, (2, 2), sl)]
    C4, c4cert = basis_from(sp4, co4)
    C3_22, _ = basis_from([npl(c) for c in spaces.spanning(3, (2, 2), sl)], co3)
    B4 = [ce.d_ce(c, 3, sl) for c in C3_22]
    vB = [co4.vec(c) for c in B4]
    dC4 = [ce.d_ce(c, 4, sl) for c in C4]
    v5 = [co5.vec(c) for c in dC4]
    k5, M5 = la.matrix_from(v5)
    if k5:
        c5 = la.certified_rank(M5)
        Z4 = [combo(C4, v) for v in c5['kernel']]
        rk5 = c5['rank']
    else:
        Z4 = list(C4)
        rk5 = 0
    keys4, MB = la.matrix_from(vB + [co4.vec(c) for c in Z4] + [co4.vec(F223_formal)])
    nB = len(vB)
    cB = la.certified_rank([row[:nB] for row in MB]) if nB else dict(rank=0, pivcols=[], pivrows=[])
    dimB4 = cB['rank']
    print("  (c) dim C^4(2,2) = %d, rank d4 = %d, dim Z^4 = %d, dim B^4 = %d, dim H^4(2,2) = %d" % (
        len(C4), rk5, len(Z4), dimB4, len(Z4) - dimB4), flush=True)
    R['H4'] = (len(C4), rk5, len(Z4), dimB4, len(Z4) - dimB4)
    # H^4 basis starting with F_223
    Bsel = [B4[j] for j in cB['pivcols']]
    chosen, cur = [], list(Bsel)
    for nm, c in [("F_223", F223_formal)] + [("z%d" % i, z) for i, z in enumerate(Z4)]:
        kk, MM = la.matrix_from([co4.vec(x) for x in cur + [c]])
        if la.certified_rank(MM)['rank'] > len(cur):
            cur.append(c)
            chosen.append((nm, c))
    hnames = [nm for nm, c in chosen]
    F_nonzero = bool(chosen) and chosen[0][0] == "F_223"
    print("  [F_223] != 0: %s;  H^4 basis: %s" % (F_nonzero, hnames))
    R['F_nonzero'] = F_nonzero
    # check [F_133 + F_223] = 0 at this N
    kk, MM = la.matrix_from([co4.vec(x) for x in Bsel] + [co4.vec(ce.lin((1, F133_formal), (1, F223_formal)))])
    if Bsel:
        cc_ = la.certified_rank([row[:-1] for row in MM])
        ex = la.in_span([row[:-1] for row in MM], [row[-1] for row in MM], cc_)[0]
    else:
        ex = all(row[-1] == 0 for row in MM)
    print("  F_133 + F_223 exact: %s" % ex)
    # coordinates

    def coords(c):
        vecs = [co4.vec(x) for x in [h for _, h in chosen] + Bsel] + [co4.vec(c)]
        kk, MM = la.matrix_from(vecs)
        A = [row[:-1] for row in MM]
        b = [row[-1] for row in MM]
        okk, x = la.in_span(A, b)
        assert okk, "bracket class not in Z^4 span (not closed?)"
        return x[:len(chosen)]
    Q = {}
    for i in range(len(Jbasis)):
        for j in range(i, len(Jbasis)):
            br = ce.bracket(Jbasis[i], 3, Jbasis[j], 3, sl=sl)
            closedN = not co5.vec(ce.d_ce(br, 4, sl))
            Q[(i, j)] = coords(br) if chosen else []
            print("  {%s,%s}: closed %s, class %s" % (Jnames[i], Jnames[j], closedN,
                  dict(zip(hnames, [str(x) for x in Q[(i, j)]]))), flush=True)
    # quadratic map and membership
    xs = sp.symbols('x0:%d' % max(1, len(Jbasis)))
    qv = [0] * len(chosen)
    for (i, j), cvec in Q.items():
        w = xs[i] * xs[j] * (sp.Rational(1, 2) if i == j else 1)
        for t in range(len(chosen)):
            qv[t] += w * sp.Rational(cvec[t].numerator, cvec[t].denominator)
    qv = [sp.factor(sp.expand(x)) for x in qv]
    print("  (b) (1/2){J,J}, J = sum x_i %s :" % Jnames, dict(zip(hnames, qv)))
    R['q'] = (Jnames, hnames, qv)
    if F_nonzero:
        y = sp.Symbol('y')
        gens = list(xs[:len(Jbasis)]) + [y]
        eqs = [e for e in qv[1:] if e != 0] + [1 - y * qv[0]]
        G = sp.groebner(eqs, *gens, order='grevlex', domain='QQ')
        member = not (len(G.exprs) == 1 and G.exprs[0] == 1)
        print("  [F_223] in image of the quadratic map over C: %s  (Groebner basis of <q_2..q_h, 1 - y q_1>: %s)" % (
            member, G.exprs if len(G.exprs) < 6 else str(len(G.exprs)) + " elements"), flush=True)
        R['member'] = member
    else:
        R['member'] = None
        print("  [F_223] = 0 at this N: [O_2] = 0, J_1 = 0 solves the order-hbar^2 equation")
    R['time'] = time.time() - t0
    print("  time %.0fs" % R['time'])
    RESULTS[name] = R
pickle.dump(RESULTS, open((__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/deform/code/algebras_%s.pkl' % ("_".join(sys.argv[1:]) or "all"), 'wb'))
print("total %.0fs" % (time.time() - t00))
