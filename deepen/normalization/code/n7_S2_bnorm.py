"""(c) side check: the raw-integral sign for S2 ((H,T) = (3,0), Omega = i dx^1 dx^2 dx^3 odd) by the same chain:
G = vs_ext * s1 * s2' * (v_3 pi^3)^{-n} K_htk,  v_3 = 8,  K_htk at (3,0) in S2 labelling,  I^{S2} = (-1)^{n(n+1)/2} K_htk
(htk CONVENTIONS Sec. 3).  s2' reorders (wedge_v Omega_v)[prod_v dxb^1_v prod_v dxb^2_v prod_v dxb^3_v] into
prod_v (Omega_v dxb^1_v dxb^2_v dxb^3_v), Omega_v counted as three odd generators.  Compared with S2 eq-B-norm:
G = (-1/(8 pi^3))^n I^{S2}."""
import sys
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/normalization/code')
import htk
from htk import conventions
from wick import perm_sign
from n3_extraction import extraction_ratio

for name, G in [("box", htk.box()), ("graph_A", htk.graph_A()), ("graph_B", htk.graph_B())]:
    vals = {extraction_ratio(G, (4, 4, 4), s)[0] for s in range(3)}
    n = len(G.cycle_matrix()); m = G.nv - 1
    vs = [v / 64 ** n for v in vals]
    s1 = (-1) ** ((3 * n) * (3 * m))
    # generators: Omega_v -> ('O', v, k) k=0,1,2 ; dxb^a_v -> ('X', a, v)
    start = [('O', v, k) for v in range(m) for k in range(3)] + [('X', a, v) for a in range(3) for v in range(m)]
    target = [g for v in range(m) for g in ([('O', v, k) for k in range(3)] + [('X', a, v) for a in range(3)])]
    pos = {g: i for i, g in enumerate(target)}
    s2 = perm_sign([pos[g] for g in start])
    mine = vs[0] * s1 * s2                       # G = mine * (8 pi^3)^{-n} K_htk
    s2kernel = (-1) ** (n * (n + 1) // 2)        # I^{S2} = s2kernel * K_htk
    s2claim = (-1) ** n * s2kernel               # S2 eq-B-norm: G = (-1)^n (8pi^3)^{-n} I^{S2}
    print("%-8s n=%d: vs_ext=%s s1=%+d s2'=%+d -> G = %+d (8pi^3)^{-n} K_htk = %+d (8pi^3)^{-n} I^S2 ; "
          "S2 eq-B-norm: %+d (8pi^3)^{-n} I^S2 ; agree: %s" % (name, n, vs, s1, s2, mine, mine * s2kernel, (-1) ** n,
                                                                 mine == s2claim))
