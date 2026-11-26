import sys, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/simple-g/code')
import algebras as al
for key in sys.argv[1:]:
    t0 = time.time()
    L = al.build(key)
    line = "%-14s dim %3d  s(kappa^-1 denominator) %d  form: %s" % (L.name, L.dim, L.s, L.form)
    if L.cartan is not None:
        C, nr, npos = al.cartan_matrix_of(L)
        line += "\n    rank %d, roots %d, Cartan matrix %s, theta = basis %d weight %s" % (len(L.cartan), nr, C.tolist(), L.theta, L.wt[L.theta])
    print(line, " (%.1fs)" % (time.time() - t0), flush=True)
