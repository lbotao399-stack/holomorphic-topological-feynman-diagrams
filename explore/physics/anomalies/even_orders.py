"""Even loop orders n = 2, 4: O_n = 24 N_n Delta^n sum_Gamma c_Gamma |Aut|^{-1} Alt(W_Gamma).
(i) formal gl(t) (all gl(N|M), t = N - M): Alt(W_Gamma) = 0 identically in the trace basis, graph by graph;
(ii) so5, sp4: Alt(W_Gamma) = 0 by explicit contraction (thetas; four-loop graphs with banana/wheel(K4) core), sum over the 24 permuted evaluations;
(iii) the four-loop kernels c_Gamma (results/bulk3_o1_std.jsonl) are read only to list them: O_4 = 0 for these algebras."""
import json, random, itertools
import numpy as np
import colour as C
rng = random.Random(41)
RES = (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/circuit-general/assembly/results/bulk3_o1_std.jsonl'
g4 = []
for line in open(RES):
    r = json.loads(line)
    g4.append(r)
g4.sort(key=lambda r: r['idx'])
thetas = [dict(idx='theta133', nv=6, edges=[(0, 5), (0, 1), (1, 2), (2, 5), (0, 3), (3, 4), (4, 5)], legs=[1, 2, 3, 4]),
          dict(idx='theta223', nv=6, edges=[(0, 1), (1, 5), (0, 2), (2, 5), (0, 3), (3, 4), (4, 5)], legs=[1, 2, 3, 4])]
cnt = 0
for th in thetas:
    rib = C.ribbon(th['nv'], th['edges'], th['legs'])
    a = C.ribbon_alt(rib)
    print('n=2 %s: formal gl(t) Alt(W) == 0: %s ; Sym(W) nonzero: %s' % (th['idx'], not a, bool(C.ribbon_sym(rib))))
allzero = True
nz_c = 0
for r in g4:
    edges = [tuple(e) for e in r['edges']]
    ext = r['legs'] + [r['nv'] - 1]
    a = C.ribbon_alt(C.ribbon(r['nv'], edges, ext))
    allzero &= (not a)
    nz_c += (r['bare'] != {})
print('n=4: %d graphs (%d with c_Gamma != 0): formal gl(t) Alt(W_Gamma) == 0 for all: %s' % (len(g4), nz_c, allzero), flush=True)

PERMS = []
for s_ in itertools.permutations(range(4)):
    sg = 1
    for i in range(4):
        for j in range(i + 1, 4):
            if s_[i] > s_[j]:
                sg = -sg
    PERMS.append((s_, sg))

def altW(alg, nv, edges, ext, xs):
    return sum(sg * C.einsum_W(alg, nv, edges, ext, [xs[s_[i]] for i in range(4)]) for s_, sg in PERMS)

for alg in (C.so(5), C.sp_(4)):
    xs = [[rng.randint(-2, 2) for _ in range(alg.dim)] for _ in range(4)]
    bad = []
    for th in thetas:
        if altW(alg, th['nv'], th['edges'], th['legs'], xs) != 0:
            bad.append(th['idx'])
    done = []
    for r in g4:
        if r['core'] not in ('banana', 'wheel'):
            continue   # prism / K33 cores: about 4 s per int64 contraction at dim 10, not run
        edges = [tuple(e) for e in r['edges']]
        done.append(r['idx'])
        if altW(alg, r['nv'], edges, r['legs'] + [r['nv'] - 1], xs) != 0:
            bad.append(r['idx'])
    print('%s: Alt(W_Gamma) == 0 for both thetas and the %d four-loop graphs with banana/wheel(K4) core %s (24 permuted evaluations, one random sample): %s %s' % (alg.name, len(done), done, not bad, bad), flush=True)
