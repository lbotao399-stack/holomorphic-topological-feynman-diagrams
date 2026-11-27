"""H = 1: hteng (rooted thetas, closed-form kernels, sequential contraction) vs the 3d N=2 engine
(explore/physics/3dN2-supercharge/ht3d.py, copied as ref_ht3d.py, labelled bitriangles, htk Taylor
tables kernels_4.pkl).  Q_0 and Q_2^hat compared on all monomials up to (deg, total derivatives).
Also: fan / tadpole parts separately.  Log: check_h1_vs_3dN2.log"""
import sys, time, itertools
from fractions import Fraction as Fr
sys.path.insert(0, ".")
import hteng, theories as T
import ref_ht3d as R
import ref_theories as RT

log = open("check_h1_vs_3dN2.log", "w")
def out(*a):
    s = " ".join(str(x) for x in a)
    print(s); log.write(s + "\n"); log.flush()

def to_ref(m):
    return tuple((n, i, a[0]) for (n, i, a) in m)

def from_ref(m):
    return tuple((n, i, (d,)) for (n, i, d) in m)

def cmp(P_mine, P_ref, A_ref, th):
    # convert mine to ref monomials (re-canonicalise in ref ordering)
    conv = {}
    for m, c in P_mine.items():
        s, mm = A_ref.canon([(n, i, a[0]) for (n, i, a) in m])
        if s:
            hteng.padd(conv, mm, s * c)
    diff = dict(conv)
    for m, c in P_ref.items():
        hteng.padd(diff, m, -c)
    return diff

cases = [
    ("U(1) charges 1,2,-3, W = phi1 phi2 phi3",
     lambda: T.u1_3d([1, 2, -3], W={((0, 0), (1, 0), (2, 0)): 1}),
     lambda: RT.u1([1, 2, -3], W={((0, 0), (1, 0), (2, 0)): 1}), 4, 2),
    ("SU(2)+adj (so3), W = 0", lambda: T.Theory(1, 3, T.F_SO3, bf=[(3, T.adj_rho(T.F_SO3, 3))]),
     lambda: RT.su2_adj(1), 4, 1),
    ("GL(2) + 1F (tadpole d^V != 0)", None, lambda: RT.gl_theory(2, 1, 0), 3, 0),
    ("3d N=4 SQED Nf=1", lambda: T.sqed_3d(1), lambda: RT.u1([0, 1, -1], W={((0, 0), (1, 0), (2, 0)): 1}), 4, 1),
]
for name, mk, mkref, maxdeg, maxd in cases:
    t0 = time.time()
    thr = mkref()
    if mk is None:
        # build GL(2)+1F in hteng from the reference data
        th = hteng.Theory(1, thr.G, thr.f, bf=[(m[0], m[1]) for m in thr.multiplets], W=thr.W)
    else:
        th = mk()
    # interaction identical?
    Iconv = {}
    for m, c in th.I.items():
        s, mm = thr.A.canon([(n, i, a[0]) for (n, i, a) in m])
        hteng.padd(Iconv, mm, s * c)
    dI = dict(Iconv)
    for m, c in thr.I.items():
        hteng.padd(dI, m, -c)
    E = hteng.Engine(th)
    ER = R.Engine(thr, R.load_kernels(4))
    pool = [L for L in th.letter_pool()]
    pool = [(n, i, (d,)) for (n, i, _) in pool for d in range(maxd + 1)]
    nm = n0 = n2 = nz2 = 0
    nfan = ntad = 0
    fanidx = [k for k, sh in enumerate(E.shapes) if "@N" in sh.name]
    tadidx = [k for k, sh in enumerate(E.shapes) if "@N" not in sh.name]
    for deg in range(1, maxdeg + 1):
        for combo in itertools.combinations_with_replacement(pool, deg):
            if sum(L[2][0] for L in combo) > maxd:
                continue
            O = th.mono(*combo)
            if not O:
                continue
            (m, c), = O.items()
            nm += 1
            q0 = E.Q0(O)
            sR, mR = thr.A.canon(to_ref(m))
            OR = {mR: c * sR}
            if cmp(q0, ER.Q0(OR), thr.A, th):
                n0 += 1
            q2 = E.Q2(O)
            nz2 += bool(q2)
            if cmp(q2, ER.Q2(OR), thr.A, th):
                n2 += 1
                if n2 <= 3:
                    out("   mismatch on", hteng.pstr(O), ":", hteng.pstr(q2), " ref:", R.pstr(ER.Q2(OR)))
            if cmp(E.Q2(O, which=fanidx), ER.Q2(OR, 'fan'), thr.A, th):
                nfan += 1
            if cmp(E.Q2(O, which=tadidx), ER.Q2(OR, 'tad'), thr.A, th):
                ntad += 1
    out("%s: I identical: %s; %d monomials (deg <= %d, derivatives <= %d): Q0 mismatches %d; "
        "Q2^hat != 0 on %d, mismatches %d (fan part %d, tadpole part %d)  [%.0fs]"
        % (name, not dI, nm, maxdeg, maxd, n0, nz2, n2, nfan, ntad, time.time() - t0))
