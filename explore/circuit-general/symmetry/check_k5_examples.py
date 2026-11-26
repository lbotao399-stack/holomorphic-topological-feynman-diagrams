"""K5 (six loops) examples at q = 7, 8: all chains found (budget 200), values compared.
Edge order (0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)."""
import engine as en

en.NUMER[0] = 3
K5 = [(i, j) for i in range(5) for j in range(i + 1, 5)]
for al, q in [((5, 5, 6, 1, 6, 4, 5, 6, 5, 5), 8), ((4, 4, 6, 7, 6, 1, 5, 5, 5, 5), 8),
              ((6, 5, 5, 4, 3, 1, 7, 7, 7, 3), 8), ((4, 2, 5, 5, 5, 6, 5, 2, 5, 3), 7)]:
    V, E, inf = en.complete(K5, al, q)
    rs = en.solve(V, E, q, budget=200, all_routes=True)
    norm = en.F(1) / en.gam(q)
    for a in al:
        norm *= en.gam(q - a)
    vals = {tuple(sorted(en.vscale(v, norm).items())) for v, ch in rs}
    print('I_%d(K5; %s) = %s   routes %d, distinct values %d' % (q, al, en.vstr(en.vscale(rs[0][0], norm)), len(rs), len(vals)))
    for v, ch in rs:
        print('    ' + ' '.join(ch))
