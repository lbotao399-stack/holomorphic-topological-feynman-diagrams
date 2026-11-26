import sys, random
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
import htk
import hf_brute as hb
import b21
nv, edges, ints = b21.banana((1, 3, 3))
nv2, edges2 = b21.replace_edge(nv, edges, 2, (1, 3, 3))
q = nv - 1
perm = list(range(nv2)); perm[q], perm[nv2 - 1] = nv2 - 1, q
E = [(perm[a], perm[b]) for a, b in edges2]
print(nv2, E)
Gm = hb.G(nv2, E); C = Gm.cycles()
rng = random.Random(1)
t = [Fr(rng.randint(1, 20), rng.randint(1, 9)) for _ in E]
P, U = hb.rho_form(Gm, t, C, list(range(len(E) - 1)))
print("rho_4 components:", len(P))
lam = [[Fr(rng.randint(-5, 5)), Fr(rng.randint(-5, 5))] for _ in range(nv2 - 1)]
y = hb.currents(Gm, t, lam)
rows = [[C[r][e] * y[e][i] for e in range(len(E) - 1)] for i in range(2) for r in range(4)]
import itertools
nz = 0
for S in itertools.combinations(range(len(E) - 1), 8):
    if hb.fdet([[row[e] for e in S] for row in rows]) != 0:
        nz += 1
print("hol 8-form nonzero components:", nz, "of", len(list(itertools.combinations(range(len(E) - 1), 8))))
c, U = hb.density(Gm, t, lam, T=1, C=C)
print("density coeff", c)
from htk import spred
G = htk.Graph(nv2, E)
red = spred.sp_reduce(G, D=3, prefer="large")
print(red.describe())
dens = red.density(lam, T=1, check=2)
print("engine num == 0:", dens.num == 0, " exps", dens.exponents())
pt = red.random_point(rng)
tt = red.tmap(pt)
c2, U2 = hb.density(Gm, tt, lam, T=1, C=C)
print("brute at engine point:", c2)
