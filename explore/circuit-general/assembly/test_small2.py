"""second moments: assembly (order 2) vs cores.kernel (exact interpolation route), C^3 three-loop K4 subdivisions."""
import sys, time
sys.path.insert(0, '../cores')
import cores
from check_c3_three_loop import K4SUB, MOM
import assembly

for name, edges, (r, z) in K4SUB[1:3]:
    t = time.time()
    res = assembly.kernel_moments(8, edges, MOM, order=2, verbose=False)
    print(name, 'order 2 assembly %.1fs, %d values' % (time.time() - t, len(res)), flush=True)
    for a, b in [((0, 0), (0, 0)), ((3, 1), (5, 2)), ((9, 0), (2, 2)), ((4, 1), (4, 2))]:
        t = time.time()
        v = cores.kernel(8, edges, MOM, [a, b])
        key = ('yy',) + tuple(sorted([a, b])[0]) + tuple(sorted([a, b])[1])
        got = res.get(key, {})
        print('   y%s y%s' % (a, b), {k: str(x) for k, x in v.items()}, 'MATCH' if v == got else 'MISMATCH %s' % got, '%.1fs' % (time.time() - t), flush=True)
