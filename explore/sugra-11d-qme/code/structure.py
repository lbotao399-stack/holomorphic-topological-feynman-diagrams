"""Part (a): decorated theta graphs of {I^12}_0 for the 11d action (S3 eq. R-action).

Flow orientation of an edge: from its (mu,nu) end to its (gamma,beta) end (S3 Section sec-assembly).
Vertex types: (a) beta d mu [1 in, 1 out], (b_k) nu^k mu mu d gamma [1 in (gamma), out: mu,mu (+nu^k)],
(c) gamma d gamma d gamma [3 in].
Rules used:
  - pairing (S3 eq-Delta): edges join mu-gamma or nu-beta;
  - D_theta (S3 eq-Ktheta, eq-insertion): D_theta is a factor of K_theta[F] and is multilinear in the
    momenta lambda_v of the 10 interior vertices, so every interior vertex needs an external letter.
Prints, for every theta_pqr, the admissible flow orientations, the number m of (c) vertices
(power of g) and the external letters; counts all 12-vertex type assignments by brute force.
"""
from itertools import product
from collections import Counter, defaultdict
import sys
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
from htk.graph import theta

PQR = [(1, 6, 6), (2, 5, 6), (3, 4, 6), (3, 5, 5), (4, 4, 5)]

# letters of each vertex type: list of (field, kind) kind in {'in','out'}; nu letters 'out' (pair with beta)
def letters(typ, k=0):
    if typ == 'a':
        return [('beta', 'in'), ('mu', 'out')]
    if typ == 'b':
        return [('gamma', 'in'), ('mu', 'out'), ('mu', 'out')] + [('nu', 'out')] * k
    if typ == 'c':
        return [('gamma', 'in')] * 3

def brute_count(G, kmax=1):
    """Brute force over vertex types (a, b_0..b_kmax, c) and over letter assignments to edge ends.
    Returns Counter[(m, external content)] of the number of admissible (types, letter-to-edge) maps
    with an external letter at every interior vertex."""
    nv = G.nv
    inc = defaultdict(list)
    for e, (a, b) in enumerate(G.edges):
        inc[a].append(e); inc[b].append(e)
    types = [('a', 0)] + [('b', k) for k in range(kmax + 1)] + [('c', 0)]
    res = Counter()
    # per vertex: list of options (type, tuple of (edge -> field contracted), external multiset)
    from itertools import permutations
    opts = []
    for v in range(nv):
        ov = []
        for t, k in types:
            L = letters(t, k)
            deg = len(inc[v])
            if len(L) < deg:
                continue
            for pick in permutations(range(len(L)), deg):
                ext = [L[i][0] for i in range(len(L)) if i not in pick]
                fields = tuple(L[i][0] for i in pick)
                ov.append((t, k, fields, tuple(sorted(ext))))
        # merge identical (t,k,fields,ext) with multiplicity
        opts.append(Counter(ov))
    partner = {'mu': 'gamma', 'gamma': 'mu', 'nu': 'beta', 'beta': 'nu'}
    # depth-first over vertices with edge consistency
    order = list(range(nv))
    def rec(i, assign, emap, mult):
        if i == nv:
            m = sum(1 for (t, k, f, ex) in assign if t == 'c')
            ext = Counter()
            for v, (t, k, f, ex) in enumerate(assign):
                ext.update(ex)
            na = sum(1 for (t, k, f, ex) in assign if t == 'a')
            res[(m, na, tuple(sorted(ext.items())))] += mult
            return
        v = order[i]
        for (t, k, fields, ex), c in opts[v].items():
            if G.degree(v) == 2 and len(ex) == 0:
                continue          # D_theta = 0: no external letter at an interior vertex
            ok = True
            new = dict(emap)
            for e, fld in zip(inc[v], fields):
                if e in new:
                    if partner[new[e]] != fld:
                        ok = False; break
                else:
                    new[e] = fld
            if ok:
                rec(i + 1, assign + [(t, k, fields, ex)], new, mult * c)
    rec(0, [], {}, 1)
    return res

def flows(G, interior):
    """Flow orientations (K = 0): interior (in,out) in {(1,1),(0,2),(2,0)}, trivalent in {(1,2),(3,0)}."""
    out = []
    for o in product((0, 1), repeat=G.ne):   # o[e]=0: flow e(0)->e(1)
        indeg = Counter(); outdeg = Counter()
        for e, (a, b) in enumerate(G.edges):
            t, h = (a, b) if o[e] == 0 else (b, a)
            outdeg[t] += 1; indeg[h] += 1
        ok = True; roles = {}
        for v in range(G.nv):
            io = (indeg[v], outdeg[v])
            if G.degree(v) == 2:
                roles[v] = {(1, 1): 'pass', (0, 2): 'source', (2, 0): 'sink'}[io]
            else:
                if io == (1, 2): roles[v] = 'triv_b'
                elif io == (3, 0): roles[v] = 'triv_c'
                else: ok = False
        if ok:
            out.append((o, roles))
    return out

if __name__ == '__main__':
    print('flow orientations with K = 0 (all (b) of type b_0):')
    for pqr in PQR:
        G, interior = theta(*pqr)
        F = flows(G, interior)
        byM = Counter()
        for o, roles in F:
            m = sum(1 for r in roles.values() if r in ('sink', 'triv_c'))
            nmu = sum(1 for r in roles.values() if r == 'pass')
            ngam = sum(1 for r in roles.values() if r in ('source', 'sink'))
            byM[(m, nmu, ngam)] += 1
        print(pqr, 'admissible flows:', len(F), ' (m, #ext mu, #ext gamma): count =', dict(sorted(byM.items())))
    print()
    print('brute force over vertex types a, b_0, b_1, c and letter-to-edge maps (theta_166, theta_445):')
    for pqr in [(1, 6, 6), (4, 4, 5)]:
        G, interior = theta(*pqr)
        R = brute_count(G, kmax=1)
        for key in sorted(R):
            print(pqr, 'm =', key[0], 'n_a =', key[1], 'external', dict(key[2]), 'maps:', R[key])
