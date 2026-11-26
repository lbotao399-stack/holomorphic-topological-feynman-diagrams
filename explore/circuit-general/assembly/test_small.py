"""assembly.kernel_moments vs cores.kernel (exact interpolation route) on C^3 three-loop K4 subdivisions."""
import sys, time
sys.path.insert(0, '../cores')
from fractions import Fraction as F
import cores
from check_c3_three_loop import K4SUB, MOM
import assembly

for name, edges, (r, z) in K4SUB[:2]:
    t = time.time()
    res = assembly.kernel_moments(8, edges, MOM, order=1, verbose=True)
    ta = time.time() - t
    ref = {k: x for k, x in (('1', r), ('z3', z)) if x}
    print(name, 'bare', res['bare'] == ref, '%.1fs' % ta, flush=True)
    nbad = 0
    for e, i in [(0, 0), (3, 1), (5, 2), (9, 0), (7, 1)]:
        t = time.time()
        v = cores.kernel(8, edges, MOM, [(e, i)])
        got = res.get(('y', e, i), {})
        print('   y[%d][%d]' % (e, i), {k: str(x) for k, x in v.items()}, 'MATCH' if v == got else 'MISMATCH %s' % got, '%.1fs' % (time.time() - t), flush=True)
