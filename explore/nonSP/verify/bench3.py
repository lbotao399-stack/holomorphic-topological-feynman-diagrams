from lattice import run
C2N3 = [(0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 4)]
lam = [[0, 1], [0, 0], [-1, -1], [1, 0]]
ex = -1/576
for h, tol in [(1.0, 1e-7), (1.0, 1e-9), (0.8, 1e-9), (1.2, 1e-10)]:
    v0, v1, b, n, dt = run(5, C2N3, [lam], h, tol)
    print("h=%.2f tol=%.0e  K0=%.13e K1=%.13e  rel0=%.2e rel1=%.2e  below=%.1e  pts=%d  %.0fs" % (h, tol, v0[0], v1[0], abs(v0[0]/ex-1), abs(v1[0]/ex-1), b, n, dt), flush=True)
