"""Independent cross-check with the run-1 exactness machinery (anomaly-structure/code/c_obstruction.run, imported
read-only): full tensors, all 9 canonical monomials of degree (2,2) for eta (no SL2 reduction), eta spanned by
trace tensors of the defining representation and f.  Structure constants and matrices from algebras.py."""
import sys, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/simple-g/code')
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/anomaly-structure/code')
import numpy as np
import algebras as al
import c_obstruction as co


class Adapter:
    def __init__(self, L):
        self.name, self.dim, self.mats = L.name, L.dim, L.mats
        self.fu, self.fl, self.kappa = L.fu, L.fl, L.kappa
        self.kinv, self.kinv_scale = L.G, L.s


for key in sys.argv[1:]:
    t0 = time.time()
    L = al.build(key)
    res, _ = co.run(Adapter(L))
    keep = {k: res[k] for k in ('d_omega133_zero', 'd_omega223_zero', 'n_eta_spanning', 'rank d(C3)',
                                'rank [d(C3) | omega]', 'exact', 'exact_133', 'exact_223')}
    print(L.name, keep, "(%.1fs)" % (time.time() - t0), flush=True)
