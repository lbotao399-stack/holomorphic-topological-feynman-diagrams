"""(S3) rho_n ^ rho_n = 0 for every n >= 1 (T >= 2 bulk kernels vanish).

A  generic symmetric M (independent dM_ab), n = 2,4,6,8: Pf(A)^Pf(A) = 0, tr (X^X)^j = 0 (X = M^-1 dM),
   A = -1/2 dM M^-1 dM; restriction to a random 2n-plane, arithmetic mod p (two primes); control Pf(A) != 0.
B  graphs: all skeleton multigraphs with n = 2,3,4 (exact Fractions, chain-sum variables);
   random cubic skeletons n = 6 (mod p) and n = 8 (mod p, 16-plane);
   chain-sum pullback: M depends on t only through chain sums (exact, admissible graphs D = 5,6, n = 2,3);
   admissible graphs D = 5, 6, n = 2 (all, exact) and random admissible subdivisions n = 4 at D = 4,5,6
   (mod p, random 8-plane).
"""
import random
import time
from fractions import Fraction as Fr

from common import Log, wedge, fadd, fscale, pfaffian, inverse
import c4graphs as cg

log = Log("c4_s3.log")
rng = random.Random(31337)
PRIMES = [2 ** 61 - 1, 1000000007]


def inv_mod(M, p):
    n = len(M)
    A = [[x % p for x in row] + [int(i == j) for j in range(n)] for i, row in enumerate(M)]
    for c in range(n):
        piv = next(r for r in range(c, n) if A[r][c] % p)
        A[c], A[piv] = A[piv], A[c]
        iv = pow(A[c][c], p - 2, p)
        A[c] = [(x * iv) % p for x in A[c]]
        for r in range(n):
            if r != c and A[r][c]:
                f = A[r][c]
                A[r] = [(x - f * y) % p for x, y in zip(A[r], A[c])]
    return [row[n:] for row in A]


def A_matrix(dM, Minv, mod):
    """A_rs = -1/2 sum_ab dM_ra Minv_ab ^ dM_bs  (dM: matrix of 1-forms)"""
    n = len(dM)
    half = (mod + 1) // 2 if mod else Fr(1, 2)
    A = [[{} for _ in range(n)] for _ in range(n)]
    for r in range(n):
        for s in range(n):
            acc = {}
            for a in range(n):
                for b in range(n):
                    c = Minv[a][b]
                    if c and dM[r][a] and dM[b][s]:
                        acc = fadd(acc, wedge(dM[r][a], dM[b][s], mod), c, mod)
            A[r][s] = fscale(acc, -half if not mod else (mod - half) % mod, mod)
    return A


def X_matrix(dM, Minv, mod):
    n = len(dM)
    return [[_lin([(Minv[r][a], dM[a][s]) for a in range(n)], mod) for s in range(n)] for r in range(n)]


def _lin(terms, mod):
    acc = {}
    for c, f in terms:
        if c and f:
            acc = fadd(acc, f, c, mod)
    return acc


def mat_form_mul(Xa, Xb, mod):
    n = len(Xa)
    return [[_sumw([(Xa[r][k], Xb[k][s]) for k in range(n)], mod) for s in range(n)] for r in range(n)]


def _sumw(pairs, mod):
    acc = {}
    for f, g in pairs:
        if f and g:
            acc = fadd(acc, wedge(f, g, mod), 1, mod)
    return acc


def trace_forms(Xm, mod):
    acc = {}
    for i in range(len(Xm)):
        acc = fadd(acc, Xm[i][i], 1, mod)
    return acc


def restrict(one_forms, nvars, k, mod):
    """substitute dt_e -> sum_j R_ej eta_j (random), k = target dimension"""
    R = [[rng.randrange(1, mod) for _ in range(k)] for _ in range(nvars)]

    def sub(f):
        out = {}
        for mask, c in f.items():
            e = mask.bit_length() - 1
            for j in range(k):
                out[1 << j] = (out.get(1 << j, 0) + c * R[e][j]) % mod
        return {a: b for a, b in out.items() if b}
    return sub


# ------------------------------------------------------------------------------------------------- A
log("A. generic symmetric M, independent dM_ab, random 2n-plane, mod p")
for n in (2, 4, 6, 8):
    for p in PRIMES[: (2 if n <= 6 else 1)]:
        dotr = n <= 6
        t0 = time.time()
        k = 2 * n
        M = [[0] * n for _ in range(n)]
        for a in range(n):
            for b in range(a, n):
                M[a][b] = M[b][a] = rng.randrange(p)
        Minv = inv_mod(M, p)
        gens = {}
        for a in range(n):
            for b in range(a, n):
                gens[(a, b)] = {1 << j: rng.randrange(1, p) for j in range(k)}
        dM = [[gens[(min(a, b), max(a, b))] for b in range(n)] for a in range(n)]
        A = A_matrix(dM, Minv, p)
        Pf = pfaffian(A, p)
        PP = wedge(Pf, Pf, p)
        trs = []
        if dotr:
            Xm = X_matrix(dM, Minv, p)
            Y = mat_form_mul(Xm, Xm, p)
            Yj = Y
            for j in range(1, n + 1):
                trs.append(len(trace_forms(Yj, p)) == 0)
                if j < n:
                    Yj = mat_form_mul(Yj, Y, p)
        trtxt = ("tr(X^X)^j = 0 for j=1..%d: %s" % (len(trs), all(trs))) if dotr else "traces not computed"
        log("  n=%d p=%d: Pf != 0: %s;  Pf^Pf = 0: %s;  %s   (%.1fs)"
            % (n, p, len(Pf) > 0, len(PP) == 0, trtxt, time.time() - t0))

# ------------------------------------------------------------------------------------------------- B
log("B. graphs")


def rho_sq_exact(nv, E):
    C = cg.cycle_matrix(nv, E)
    assert cg.check_cycles(nv, E, C)
    n = len(C)
    t = [Fr(rng.randint(1, 30), rng.randint(1, 7)) for _ in E]
    M = [[sum(C[a][e] * C[b][e] * t[e] for e in range(len(E))) for b in range(n)] for a in range(n)]
    Minv = inverse(M)
    dM = [[{1 << e: C[a][e] * C[b][e] for e in range(len(E)) if C[a][e] * C[b][e]} for b in range(n)]
          for a in range(n)]
    A = A_matrix(dM, Minv, None)
    Pf = pfaffian(A)
    return len(Pf) > 0, len(wedge(Pf, Pf)) == 0


def rho_sq_mod(nv, E, p, k=None):
    C = cg.cycle_matrix(nv, E)
    assert cg.check_cycles(nv, E, C)
    n = len(C)
    t = [rng.randrange(1, p) for _ in E]
    M = [[sum(C[a][e] * C[b][e] * t[e] for e in range(len(E))) % p for b in range(n)] for a in range(n)]
    Minv = inv_mod(M, p)
    if k is None:
        sub = lambda f: f
    else:
        sub = restrict(None, len(E), k, p)
    dM = [[sub({1 << e: (C[a][e] * C[b][e]) % p for e in range(len(E)) if C[a][e] * C[b][e]})
           for b in range(n)] for a in range(n)]
    A = A_matrix(dM, Minv, p)
    Pf = pfaffian(A, p)
    return len(Pf) > 0, len(wedge(Pf, Pf, p)) == 0


for n in (2, 3, 4):
    t0 = time.time()
    sks = cg.skeletons(n)
    res = [rho_sq_exact(s, E) for s, E in sks]
    nontriv = sum(1 for s, E in sks if len(E) >= 2 * n + 1)
    log("  skeletons n=%d: %d multigraphs (%d with m >= 2n+1);  rho_n != 0 on %d;  rho^rho = 0 on all: %s   (%.1fs)"
        % (n, len(sks), nontriv, sum(r[0] for r in res), all(r[1] for r in res), time.time() - t0))

for n, k, cnt in ((6, None, 6), (8, 16, 1)):
    t0 = time.time()
    res = []
    for _ in range(cnt):
        s, E = cg.random_cubic_skeleton(n, rng)
        res.append(rho_sq_mod(s, E, PRIMES[0], k))
    log("  random cubic skeletons n=%d (m=%d) x%d, mod p%s: rho != 0: %s; rho^rho = 0: %s   (%.1fs)"
        % (n, 3 * (n - 1), cnt, "" if k is None else ", 16-plane", all(r[0] for r in res), all(r[1] for r in res),
           time.time() - t0))

# chain-sum pullback: within every chain all cycle columns agree up to sign => M(t) = M_skeleton(R_c = sum_chain t)
log("  chain-sum pullback (C_ae C_be constant along every chain), admissible graphs:")


def chains(nv, E):
    deg = [0] * nv
    for a, b in E:
        deg[a] += 1
        deg[b] += 1
    inc = {v: [] for v in range(nv)}
    for e, (a, b) in enumerate(E):
        inc[a].append(e)
        inc[b].append(e)
    seen, out = set(), []
    for e0 in range(len(E)):
        if e0 in seen:
            continue
        ch, stack = [], [e0]
        while stack:
            e = stack.pop()
            if e in seen:
                continue
            seen.add(e)
            ch.append(e)
            for v in E[e]:
                if deg[v] == 2:
                    stack.extend(f for f in inc[v] if f not in seen)
        out.append(ch)
    return out


for D, n in ((5, 2), (6, 2), (5, 3), (6, 3)):
    t0 = time.time()
    gs = cg.admissible(D, n)
    ok = True
    nch = set()
    for nv, E in gs:
        C = cg.cycle_matrix(nv, E)
        chs = chains(nv, E)
        nch.add(len(chs))
        for ch in chs:
            prods = {tuple(C[a][e] * C[b][e] for a in range(n) for b in range(n)) for e in ch}
            ok &= (len(prods) == 1)
    log("    D=%d n=%d: %d graphs, chain counts %s; C_ae C_be constant on chains: %s   (%.1fs)"
        % (D, n, len(gs), sorted(nch), ok, time.time() - t0))

for D in (5, 6):
    gs = cg.admissible(D, 2)
    res = [rho_sq_exact(nv, E) for nv, E in gs]
    log("  admissible D=%d n=2: %d graphs (thetas), exact: rho_2 != 0: %s, rho_2^rho_2 = 0: %s"
        % (D, len(gs), all(r[0] for r in res), all(r[1] for r in res)))

for D in (4, 5, 6):
    t0 = time.time()
    N = 4 * D + 1
    sks = [x for x in cg.skeletons(4) if len(x[1]) >= 9]
    picked = []
    tries = 0
    while len(picked) < 8 and tries < 200000:
        tries += 1
        s, sk = rng.choice(sks)
        m = len(sk)
        cuts = sorted(rng.sample(range(1, N), m - 1))
        lens = tuple(b - a for a, b in zip([0] + cuts, cuts + [N]))
        if not cg.local_ok(D, s, sk, lens):
            continue
        nv, E = cg.subdivide(s, sk, lens)
        if len(set(tuple(sorted(e)) for e in E)) < len(E):
            continue
        picked.append((nv, E))
    res = [rho_sq_mod(nv, E, PRIMES[0], 8) for nv, E in picked]
    log("  admissible D=%d n=4 (|I|=%d): %d random subdivisions of 6-vertex skeletons, mod p, 8-plane:"
        " rho_4 != 0: %s, rho_4^rho_4 = 0: %s   (%.1fs)"
        % (D, N, len(picked), all(r[0] for r in res), all(r[1] for r in res), time.time() - t0))
log.close()
