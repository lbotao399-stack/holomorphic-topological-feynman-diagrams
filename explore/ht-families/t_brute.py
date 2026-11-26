import sys, random
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
import htk
from htk import brute
import hf_brute as hb
edges = [(0,1),(0,2),(0,4),(0,5),(1,2),(1,3),(1,4),(3,4),(4,5)]
rng = random.Random(1)
Ge = htk.Graph(6, edges); Gm = hb.G(6, edges)
for trial in range(3):
    t = [Fr(rng.randint(1, 20), rng.randint(1, 9)) for _ in edges]
    lam = [[Fr(rng.randint(-9, 9))] for _ in range(5)]
    c, U = hb.density(Gm, t, lam, T=1)
    De = brute.density(Ge, t, lam, T=1)
    print(trial, "mine", c, U, "engine", De, "same:", (c > 0) == (De.coeff > 0) and c**2 / U == De.square())
