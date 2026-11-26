"""Exact bare kernels K(lambda;0) of the bulk3(3,3) graphs (hCS on C^3, three loops, four legs).
Sources: legs (three degree-2 vertices) = e_1,e_2,e_3, fourth degree-2 vertex = root = -(e1+e2+e3) (std_sources),
plus a random integer check K(lambda) = c det(lambda_legs)^3 at a second source set.
Output: kernels_n3.json"""
import json, os, sys, random
from fractions import Fraction as Fr
A = (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/circuit-general/assembly'
sys.path.insert(0, A)
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/circuit-general/atlas')
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
os.environ.setdefault('NUMBA_CACHE_DIR', '/private/tmp/claude-501/-Users-libotao-Desktop/07c118b5-2480-4b3c-899c-a49c36ac4d00/scratchpad/numba_anom')
import networkx as nx
from atlas import set_bulk3
from htk.graph import from_networkx
import htk
import assembly
from run_batch import std_sources

def det3(a, b, c):
    return (a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0]) + a[2]*(b[0]*c[1]-b[1]*c[0]))

n = int(sys.argv[1]) if len(sys.argv) > 1 else 3
rng = random.Random(11)
out = []
for idx, it in enumerate(set_bulk3(3, n)):
    G = it[0]
    S0 = [v for v in G if G.degree(v) == 2]
    O = S0[-1]
    Gh, lab = from_networkx(G, root=O)
    legs, full = std_sources(Gh, lab, S0, O, 3)
    res = assembly.kernel_moments(Gh.nv, Gh.edges, full, order=0)
    bare = {k: str(v) for k, v in res['bare'].items()}
    # second source set
    vecs = [[rng.randint(-4, 4) for _ in range(3)] for _ in range(3)]
    full2 = [[0]*3 for _ in range(Gh.nv)]
    for a, v in enumerate(legs):
        full2[v] = vecs[a]
    full2[Gh.nv-1] = [-sum(r[i] for r in full2[:-1]) for i in range(3)]
    res2 = assembly.kernel_moments(Gh.nv, Gh.edges, full2, order=0)
    d = det3(*vecs)
    cov = all(Fr(res2['bare'].get(k, 0)) == Fr(res['bare'].get(k, 0)) * d**n for k in set(res['bare']) | set(res2['bare']))
    sp = htk.catalog.is_sp(nx.Graph(G))
    rec = dict(idx=idx, nv=Gh.nv, edges=[list(e) for e in Gh.edges], legs=legs, root=Gh.nv-1, sp=bool(sp),
               bare=bare, cov_check=bool(cov), det2=d)
    if sp:
        lam = {v: full[v] for v in range(Gh.nv-1)}
        kx = htk.kernel_exact(Gh, lam)
        rec['htk_spred'] = str(kx)
        rec['htk_match'] = (Fr(str(kx)) == Fr(res['bare'].get('1', 0))) and all(k == '1' for k in res['bare'])
    print(json.dumps(rec), flush=True)
    out.append(rec)
json.dump(out, open('kernels_n%d.json' % n, 'w'), indent=1)
