"""B3: U(N) = SU(N) x U(1) in the tree-level twisted complex, sector by sector, and SU(2) per-level
totals against the census (data/N4_SU2_levels.csv) and U(2) totals against Chang-Lin Table 3.

H^{U(N)}(lam, Y) = sum_{lam'+lam''=lam, Y'+Y''=Y} H^{SU(N)}(lam', Y') * dim Sym(V_1)(lam'', Y'')
usage: python3 b3_U_vs_SU.py N Jmax
"""
import sys
from itertools import permutations
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
from holography import bg_words as bw, uN_q0 as uq, census_io as cio

N = int(sys.argv[1])
J = int(sys.argv[2])
alg, keys = bw.n4_letters(J)


def sectors(j):
    out = []
    for n in range(j // 3 + 1):
        r = j - 3 * n
        if r % 2:
            continue
        Th = r // 2
        for npl in range(n + 1):
            for t1 in range(Th + 1):
                for t2 in range(Th - t1 + 1):
                    out.append((t1, t2, Th - t1 - t2, npl, n - npl))
    return out


def rep(lam):
    t = tuple(sorted(lam[:3], reverse=True))
    n = tuple(sorted(lam[3:], reverse=True))
    return t + n


cache = {}


def H(lam, special, NN):
    key = (rep(lam), special, NN)
    if key not in cache:
        cache[key] = {Y: v[2] for Y, v in uq.cohomology_sector(alg, rep(lam), NN, special=special).items() if v[2]}
    return cache[key]


# U(1) Fock space of letters (free): dim Sym(V_1)(lam, Y) = H^{U(1)}
allsec = [lam for j in range(0, J + 1) for lam in (sectors(j) if j else [(0, 0, 0, 0, 0)])]
fock = {}
for lam in allsec:
    fock[lam] = {0: 1} if lam == (0, 0, 0, 0, 0) else H(lam, False, 1)

cen = cio.su2_levels()
CL = cio.chang_lin_tables()
bad = []
for j in range(2, J + 1):
    BF_su = [0, 0]
    tot_u = 0
    for lam in sectors(j):
        Th = sum(lam[:3])
        hu = H(lam, False, N)
        hs = H(lam, True, N)
        for Y, h in hs.items():
            BF_su[(Y + Th) % 2] += h
        tot_u += sum(hu.values())
        # convolution
        conv = {}
        for lam1 in allsec:
            if all(a <= b for a, b in zip(lam1, lam)):
                lam2 = tuple(b - a for a, b in zip(lam1, lam))
                hs1 = {0: 1} if lam1 == (0, 0, 0, 0, 0) else H(lam1, True, N)
                f2 = fock.get(lam2, {})
                for Y1, a in hs1.items():
                    for Y2, b in f2.items():
                        conv[Y1 + Y2] = conv.get(Y1 + Y2, 0) + a * b
        conv = {Y: v for Y, v in conv.items() if v}
        if conv != hu:
            bad.append((lam, hu, conv))
    line = f"level {j:2d}: SU({N}) (B,F) = {tuple(BF_su)}"
    if N == 2 and j in cen:
        line += f"  census (B,F) = ({cen[j]['bosonic']},{cen[j]['fermionic']})"
    if N in CL and j < len(CL[N]['d_u']):
        line += f"  U({N}) total = {tot_u}  Chang-Lin d_u = {CL[N]['d_u'][j]}  d_su = {CL[N]['d_su'][j]}"
    print(line, flush=True)
print("sectors where H^{U(N)} != H^{SU(N)} * Fock(U(1)):", len(bad))
for b in bad[:10]:
    print(b)
