import sys, time
import valg
P1, P2 = 1000000007, 998244353
for name in sys.argv[1:]:
    t0 = time.time()
    A = valg.build(name)
    ch = A.checks()
    extra = {}
    if name == 'G2':
        extra = {'octonion alternative+normed': A.octonion_ok, 'dim Der(O)': A.n_der}
    if name == 'F4':
        extra = {'Clifford relations + symmetric gammas': A.clifford_ok, 'rho representation': A.rho_ok}
    nul = (valg.invariant_forms_nullity(A, P1), valg.invariant_forms_nullity(A, P2))
    rk = valg.rank_of_algebra(A)
    print("%-5s dim %d ; %s ; %s ; %s ; dim inv. sym. forms mod p1,p2 = %s ; dim ker ad_X (random X) = %d ; s = %d ; (%.1fs)"
          % (A.name, A.n, A.info, ch, extra, nul, rk, A.s, time.time() - t0), flush=True)
