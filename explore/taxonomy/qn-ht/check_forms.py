"""5d N=1 with U = T^*R: ESW Thm 9.2 (U//g, u-u propagator) vs Thm 9.3 (BF R/g, phi-psi propagator).
Map u^i -> phi^i, u^{dim+i} -> sigma psi_i (sigma = +-1); compare Q_0 and Q_2^hat on all monomials of
degree <= 3, total derivatives <= maxd.  R = fundamental of sl2 and U(1) charge 2.  Log: check_forms.log"""
import sys, itertools, time
sys.path.insert(0, ".")
import hteng, theories as T

log = open("check_forms.log", "w")
def out(*a):
    s = " ".join(str(x) for x in a)
    print(s); log.write(s + "\n"); log.flush()

def conv(P, th_to, dim, sigma):
    o = {}
    for m, c in P.items():
        lets, sg = [], c
        for (n, i, a) in m:
            if n == 'u0':
                if i < dim:
                    lets.append(('p0', i, a))
                else:
                    lets.append(('s0', i - dim, a)); sg *= sigma
            else:
                lets.append((n, i, a))
        s, mm = th_to.canon(lets)
        if s:
            hteng.padd(o, mm, s * sg)
    return o

for gname, G, F, rho, dim in [("SU(2) fund", 3, T.F_SL2, T.FUND_SL2, 2), ("U(1) q=2", 1, {}, [{(0, 0): 2}], 1)]:
    th_bf = hteng.Theory(2, G, F, bf=[(dim, rho)])
    gauge = 'sl2' if G == 3 else 'u1'
    th_u = T.sym_theory(2, gauge, [T.full_hyper_sym(rho, dim)], "U//g")
    Ebf, Eu = hteng.Engine(th_bf), hteng.Engine(th_u)
    for sigma in (1, -1):
        t0 = time.time()
        pool = [(n, i, a) for (n, i, _) in th_u.letter_pool() for a in itertools.product(range(2), repeat=2) if sum(a) <= 1]
        n = bad0 = bad2 = nz = 0
        for deg in (1, 2, 3):
            for combo in itertools.combinations_with_replacement(pool, deg):
                if sum(sum(L[2]) for L in combo) > 1:
                    continue
                O = th_u.mono(*combo)
                if not O:
                    continue
                n += 1
                Ob = conv(O, th_bf, dim, sigma)
                d0 = hteng.psum(conv(Eu.Q0(O), th_bf, dim, sigma), {m: -c for m, c in Ebf.Q0(Ob).items()})
                q2 = Eu.Q2(O)
                nz += bool(q2)
                d2 = hteng.psum(conv(q2, th_bf, dim, sigma), {m: -c for m, c in Ebf.Q2(Ob).items()})
                bad0 += bool(d0); bad2 += bool(d2)
        out("%s, sigma = %+d: %d monomials (deg <= 3, derivatives <= 1): Q0 mismatches %d; Q2^hat != 0 on %d, mismatches %d  [%.0fs]"
            % (gname, sigma, n, bad0, nz, bad2, time.time() - t0))
