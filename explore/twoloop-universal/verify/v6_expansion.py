"""C9: Laplace-expansion signs eps0, eps1 (as stated), term counts, and the data tables of results.json
(evaluated as data against own closed values with random rational sources)."""
import sys, json, random, math, itertools, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/twoloop-universal/verify")
sys.dont_write_bytecode = True
from vcore import *

rng = random.Random(1357)
DIR = (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/twoloop-universal/"
t0 = time.time()


def br(vecs):
    H = len(vecs)
    return det([[Fr(vecs[c][i]) for c in range(H)] for i in range(H)])


def eps1(p, H, Bp):
    return sg(H * (H + 1) // 2 + p * (p - 1) // 2 + (H - p + 1) * (p - 1) + sum(Bp))


def eps0(p, H, Bp):
    return sg(H * (H + 1) // 2 + (p + 1) * (p + 2) // 2 + (H - p) * (p + 1) + sum(Bp))


def Dtheta_expanded(th, lam, H):
    p, q = th.n["A"], th.n["B"]
    LA = [lam[u] for u in th.interior["A"]]
    LB = [lam[u] for u in th.interior["B"]]
    LC = [lam[u] for u in th.interior["C"]]
    k = H - p + 1
    if k < 0 or k > q - 1:
        return Fr(0), 0
    tot, n = Fr(0), 0
    for Bp in itertools.combinations(range(1, q), k):
        Bpp = [b for b in range(1, q) if b not in Bp]
        tot += eps1(p, H, Bp) * br(LA + [LB[b - 1] for b in Bp]) * br([LB[b - 1] for b in Bpp] + LC)
        n += 1
    return tot, n


def Q_expanded(th, lam, H, v):
    p, q = th.n["A"], th.n["B"]
    LA = [lam[u] for u in th.interior["A"]]
    LB = [lam[u] for u in th.interior["B"]]
    LC = [lam[u] for u in th.interior["C"]]
    k = H - p
    if k < 0 or k > q - 1:
        return Fr(0), 0
    tot, n = Fr(0), 0
    for Bp in itertools.combinations(range(1, q), k):
        Bpp = [b for b in range(1, q) if b not in Bp]
        tot += eps0(p, H, Bp) * br([v] + LA + [LB[b - 1] for b in Bp]) * br([v] + [LB[b - 1] for b in Bpp] + LC)
        n += 1
    return tot, n


# ---- eps formulas on all tight graphs D <= 7 and all orderings of (p, q, r)
ok_eps, cnt = True, 0
for D in range(1, 8):
    for pqr in theta_list(D):
        for perm in sorted(set(itertools.permutations(pqr))):
            th = theta_std(*perm)
            for T in (0, 1):
                H = D - T
                if H < 1 and T == 0:
                    continue
                lam = rand_sources(th.nv, th.Np, H, rng)
                if T == 1:
                    ex, _ = Dtheta_expanded(th, lam, H)
                    ok_eps &= (ex == D_theta(th, lam, H))
                else:
                    v = [Fr(rng.randint(-9, 9), 5) for _ in range(H)]
                    ex, _ = Q_expanded(th, lam, H, v)
                    ok_eps &= (ex == Qform(th, lam, H, v))
                cnt += 1
print("eps1 / eps0 Laplace formulas = direct 2H x 2H determinants: %s (%d graph/T cases, D <= 7, all path orders)" % (ok_eps, cnt))

# ---- S3 signs: D_355 (-1)^d, D_445 (-1)^{d+1}  (B'' = {d})
s355 = all(eps1(3, 5, [b for b in range(1, 5) if b != d]) == sg(d) for d in range(1, 5))
s445 = all(eps1(4, 5, [b for b in range(1, 4) if b != d]) == sg(d + 1) for d in range(1, 4))
print("eps1 reproduces S3 signs (-1)^d (D_355): %s, (-1)^{d+1} (D_445): %s" % (s355, s445))


# ---- term counts: canonical products of two brackets (ordered pair top/bottom), names sorted
def names(th):
    nm = {th.N: "N", th.Np: "Np"}
    for P in "ABC":
        for k, u in enumerate(th.interior[P], start=1):
            nm[u] = P + str(k)
    return nm


def canon(lst):
    if len(set(lst)) < len(lst):
        return 0, None
    order = sorted(range(len(lst)), key=lambda i: lst[i])
    return perm_sign(order), tuple(lst[i] for i in order)


def expand_terms(th, H, T):
    nm = names(th)
    p, q = th.n["A"], th.n["B"]
    A = [nm[u] for u in th.interior["A"]]
    B = [nm[u] for u in th.interior["B"]]
    C = [nm[u] for u in th.interior["C"]]
    terms = {}
    if T == 1:
        k = H - p + 1
        for Bp in itertools.combinations(range(1, q), k):
            Bpp = [b for b in range(1, q) if b not in Bp]
            s1, t1 = canon(A + [B[b - 1] for b in Bp])
            s2, t2 = canon([B[b - 1] for b in Bpp] + C)
            if s1 and s2:
                terms[(t1, t2)] = terms.get((t1, t2), Fr(0)) + eps1(p, H, Bp) * s1 * s2
    else:
        mu = {"N": Fr(1)}
        for P in "ABC":
            for k_, u in enumerate(th.interior[P], start=1):
                mu[nm[u]] = Fr(th.n[P] - k_, th.n[P])
        k = H - p
        for Bp in itertools.combinations(range(1, q), k):
            Bpp = [b for b in range(1, q) if b not in Bp]
            for a in mu:
                for b in mu:
                    s1, t1 = canon([a] + A + [B[x - 1] for x in Bp])
                    s2, t2 = canon([b] + [B[x - 1] for x in Bpp] + C)
                    if s1 and s2:
                        terms[(t1, t2)] = terms.get((t1, t2), Fr(0)) + eps0(p, H, Bp) * s1 * s2 * mu[a] * mu[b]
    return {k: v for k, v in terms.items() if v != 0}


claimed_T0 = [1, 4, 9, 9, 16, 16, 30, 25, 25, 25, 69, 36, 36, 36, 132, 102]
claimed_T1 = {(3, 3, 3): 2, (3, 4, 4): 3, (3, 5, 5): 4, (4, 4, 5): 3}
c0s, c1s = [], {}
for T in (0, 1):
    for D in range(1, 7):
        H = D - T
        if T == 0 and H < 1:
            continue
        for pqr in theta_list(D):
            n = len(expand_terms(theta_std(*pqr), H, T))
            if T == 0:
                c0s.append(n)
            else:
                c1s[pqr] = n
print("own term counts T=0:", c0s, " claim:", claimed_T0, " equal:", c0s == claimed_T0)
print("own term counts T=1:", c1s, " claim (1 unless listed):", claimed_T1,
      " equal:", all(c1s[k] == claimed_T1.get(k, 1) for k in c1s))

# ---- results.json tables as data, against own closed values
res = json.load(open(DIR + "results.json"))
okt, nt = True, 0
for tab in res["tables"]:
    T, H, pqr = tab["T"], tab["H"], tuple(tab["pqr"])
    th = theta_std(*pqr)
    nm = names(th)
    lab = {v: k for k, v in nm.items()}
    assert {k: v for k, v in tab["vertex_labels_htk_theta"].items()} == {nm[v]: v for v in nm}
    q = pqr[1]
    D = H + T
    own_terms = expand_terms(th, H, T)
    file_terms = {(tuple(tp), tuple(bt)): Fr(cf) for cf, tp, bt in tab["terms"]}
    okt &= (own_terms == file_terms)
    for _ in range(2):
        lam = rand_sources(th.nv, th.Np, H, rng)
        fac = Fr(1)
        for n_ in pqr:
            fac /= math.factorial(n_ - 1)
        if T == 1:
            c = fac * ghalf(H + 1 - pqr[0]) * ghalf(H + 1 - pqr[1]) * ghalf(H + 1 - pqr[2]) / ghalf(H + 1) / 2
            own = sg(H * (H + 1) // 2) * sg((q + 1) * D) * c * D_theta(th, lam, H)          # coefficient of pi
            ev = sum(Fr(cf.replace("*pi", "")) * br([lam[lab[x]] for x in tp]) * br([lam[lab[x]] for x in bt])
                     for cf, tp, bt in tab["K_htk_terms"])
        else:
            c = fac * Fr(math.prod(math.factorial(H - n_) for n_ in pqr), math.factorial(H + 1))
            mu = {th.N: Fr(1), th.Np: Fr(0)}
            for P in "ABC":
                for k, u in enumerate(th.interior[P], start=1):
                    mu[u] = Fr(th.n[P] - k, th.n[P])
            Lbar = [sum(mu[v] * Fr(lam[v][i]) for v in mu) for i in range(H)]
            own = sg(H * (H + 1) // 2 + 1) * sg((q + 1) * D) * c * Qform(th, lam, H, Lbar)
            ev = sum(Fr(cf) * br([lam[lab[x]] for x in tp]) * br([lam[lab[x]] for x in bt])
                     for cf, tp, bt in tab["K_htk_terms"])
            okt &= all(Fr(tab["Lbar_weights"][nm[v]]) == mu[v] for v in mu if v != th.Np)
        okt &= (ev == own) and (Fr(tab["c"]) == c) and (tab["sign_htk"] == (sg(H * (H + 1) // 2 + (1 - T)) * sg((q + 1) * D)))
    nt += 1
print("results.json tables (%d): terms = own expansion, K_htk_terms = own K_htk at random sources, c, sign_htk, Lbar weights: %s"
      % (nt, okt))
print("time %.1f s" % (time.time() - t0))
