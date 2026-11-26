import sys, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/simple-g/code')
import algebras as al, invariants as inv
SPECS = [(2, 'S'), (2, 'L'), (3, 'S'), (3, 'L'), (3, 'T'), (4, 'S'), (4, 'L'), (4, 'T')]
mode = sys.argv[1]
for key in sys.argv[2:]:
    t0 = time.time()
    L = al.build(key)
    if mode == 'weyl':
        dims, nW = inv.invariant_dims(L, SPECS)
        print("%-14s |W| = %4d  " % (L.name, nW) + "  ".join("dim(%s^%d g)^g = %d" % ({'S': 'S', 'L': 'Lambda', 'T': 'T'}[kd], k, dims[(k, kd)])
              for (k, kd) in SPECS) + "   (%.1fs)" % (time.time() - t0), flush=True)
    else:
        for (k, kd) in [(3, 'S'), (3, 'L')] + ([(2, 'L')] if mode == 'lin' else []):
            nul, nu, ne = inv.invariant_tensors_linalg(L, k, kd)
            print("%-14s linear algebra: %s^%d: unknowns %d, equations %d, nullity mod p1, p2 = %s   (%.1fs)"
                  % (L.name, {'S': 'S', 'L': 'Lambda'}[kd], k, nu, ne, nul, time.time() - t0), flush=True)
