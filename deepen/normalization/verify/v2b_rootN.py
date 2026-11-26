"""V2b: pointwise identity (V2 code) for the msym-odd theta_rootN labelling at (2,1)."""
import sys
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/normalization/verify')
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
import htk
from v2_pointwise import run


def rootN(p, q, r):
    nv = p + q + r - 1
    N, Np = nv - 1, nv - 2
    edges, nxt = [], 0
    for n_ in (p, q, r):
        prev = N
        for j in range(n_):
            cur = Np if j == n_ - 1 else nxt
            if j < n_ - 1:
                nxt += 1
            edges.append((prev, cur))
            prev = cur
    return htk.Graph(nv, edges, 'rootN')


print("V2b (2,1) theta_rootN labelling (msym-odd):")
for pqr in [(1, 3, 3), (2, 2, 3)]:
    run(rootN(*pqr), '(2,1)', 'theta_rootN%s' % (pqr,), npts=2)
