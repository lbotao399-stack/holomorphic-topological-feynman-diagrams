"""t3_banana.py -- H = 0, T = 1 (topological quantum mechanics): banana graphs, N = n+1 parallel edges.

 (a) parallel move with source kappa (vertex 0, -kappa at the root), exact at random rational t:
       <prod psi> = <(-w_*)^n> t_red^n prod_i r_i^{-1} dr_2 ^ ... ^ dr_N,
       <w_*^n> = E[(kappa + xi/sqrt(t_red))^n],  xi ~ N(0, 1/2),   U = t_red^n / prod r_i,   P = t_red kappa^2,
     hence  T[kappa] = (-1)^n e^{t_red kappa^2} E[(sqrt(t_red) kappa + xi)^n] prod_i r_i^{-1/2} dr_2...dr_N.
 (b) htk density at kappa = 0:  D(t) = eps_N E[xi^n] t_red^{N/2} prod_e t_e^{-3/2}  (sign eps_N and square, exact).
 (c) K_N = int_{sum t = 1} D = eps_N E[xi^n] Gamma(1/2)^N / Gamma(N/2)  ( = eps_N pi^{n/2}, n even; 0, n odd );
     GKW: I_N = (-1)^n (4 pi)^{-n/2} K_N ( = 2^{1-N} for N odd, Higher.tex (3.28) ).
 (d) numerical oracle: lattice rule (N = 3) and importance Monte Carlo (N = 5) in log Schwinger coordinates
     of the brute-force density (Pfaffian, floating point).
 (e) ordered product (Higher.tex (4.24)): I*_N = I_N / 2 + I_{N+1} = 2^{-N}; Moyal coefficients I*_N / N!.
 (f) UV boundary of [eps, oo)^N: t_red <= eps on every face, so the kappa-dependence of (a) drops as eps -> 0.
"""
import random, math, time
from fractions import Fraction as Fr
import tb
from tb import TopoData, wedge, fadd, fscale, pullback, bits
from wick import Dual
import graphs as gr

rng = random.Random(2718)


def dfact(m):
    out = 1
    while m > 1:
        out *= m
        m -= 2
    return out


def E_shift(a2, n, kappa, tred):
    """E[(kappa + xi/sqrt(tred))^n], xi ~ N(0,1/2) -> sum_j C(n,2j) kappa^{n-2j} (2j-1)!! / (2 tred)^j"""
    tot = Fr(0)
    for j in range(0, n // 2 + 1):
        tot += Fr(math.comb(n, 2 * j)) * Fr(kappa) ** (n - 2 * j) * dfact(2 * j - 1) / (2 * tred) ** j
    return tot


def htk_density(Omega, t, N):
    """D with Omega = D iota_E(dt_0...dt_{N-1}); e0 = 0"""
    mask = ((1 << N) - 1) ^ 1
    return Omega.get(mask, Fr(0)) / Fr(t[0])     # (-1)^0 / t_0


def main():
    t0 = time.time()
    print(__doc__)
    Kval = {}
    for N in range(2, 8):
        n = N - 1
        G = gr.banana(N)
        ok_a, ok_b, signs = True, True, set()
        for rep in range(3):
            t = tb.rand_t(N, rng)
            kap = Fr(rng.randint(-7, 7), rng.randint(1, 4))
            D = TopoData(G, t)
            Phi, U, P = D.topo_form([kap])
            T = [Dual.var(t[e], e) for e in range(N)]
            tred = 1 / sum((1 / x for x in T), Dual(0))
            r = [tred / x for x in T]
            vol = {0: Fr(1)}
            for i in range(1, N):
                vol = wedge(vol, r[i].form())
            prodr = Fr(1)
            for x in r:
                prodr *= x.v
            pred = fscale(vol, (-1) ** n * E_shift(None, n, kap, tred.v) * tred.v ** n / prodr)
            ok_a &= (Phi == pred) and (U == tred.v ** n / prodr) and (P == tred.v * kap ** 2)
            # kappa = 0 density
            Phi0, U0, _ = D.topo_form([Fr(0)])
            d0 = htk_density(Phi0, t, N)                     # D(t) = d0 * U^{-1/2}
            if n % 2 == 0:
                c = Fr(dfact(n - 1), 2 ** (n // 2))
                prodt = Fr(1)
                for x in t:
                    prodt *= x
                # D^2 = d0^2 / U  ==  c^2 tred^N prod t^-3
                ok_b &= (d0 ** 2 / U0 == c ** 2 * tred.v ** N / prodt ** 3)
                signs.add((d0 > 0) - (d0 < 0))
            else:
                ok_b &= (d0 == 0)
        eps = signs.pop() if len(signs) == 1 else None
        # exact kernel: eps * E[xi^n] * Gamma(1/2)^N / Gamma(N/2); for n even = eps * pi^{n/2}
        if n % 2 == 0:
            Kq = (eps, n // 2)   # K = eps * pi^{n/2}
            I = Fr(eps * (-1) ** n, 4 ** (n // 2))     # (-1)^n (4 pi)^{-n/2} K
        else:
            Kq, I = (0, 0), Fr(0)
        Kval[N] = I
        print("N=%d (n=%d): (a) form with source = parallel-move formula: %s; (b) density: %s, sign eps_N = %s; "
              "K_N = %s; I_N = %s" % (N, n, ok_a, ok_b, eps, ("%d pi^%d" % Kq) if n % 2 == 0 else "0", I))
    Kval[1] = Fr(1)   # single edge: alpha = 1 (Higher.tex (3.55))
    print("\n(e) ordered product coefficients I*_N = I_N/2 + I_{N+1}  (I_1 = 1):")
    for N in range(1, 7):
        Is = Kval[N] / 2 + Kval[N + 1]
        print("   N=%d: I*_N = %s, 2^-N = %s, equal: %s;  Moyal coefficient I*_N/N! = %s = (1/2)^N/N!: %s"
              % (N, Is, Fr(1, 2 ** N), Is == Fr(1, 2 ** N), Is / math.factorial(N), Is / math.factorial(N) == Fr(1, 2 ** N * math.factorial(N))))
    print("   commutator coefficients I_N/N!: " + ", ".join("N=%d: %s" % (N, Kval[N] / math.factorial(N)) for N in range(1, 8)))
    # (d) oracle: htk lattice rule (CONVENTIONS sec 6): K = int_{R^{N-1}} S(e^s) d^{N-1}s, S = D prod_e t_e,
    #     s_{N-1} = 0; S from the brute-force Pfaffian in floating point.
    print("\n(d) oracle K = int S(e^s) d^{N-1}s, S = D prod t (brute-force Pfaffian, floats):")

    def S_of(svec, N):
        t = [Fr(math.exp(x)) for x in svec] + [Fr(1)]
        D = TopoData(gr.banana(N), t)
        Phi0, U0, _ = D.topo_form([Fr(0)])
        prodt = 1.0
        for x in t:
            prodt *= float(x)
        return float(htk_density(Phi0, t, N)) / math.sqrt(float(U0)) * prodt
    # N = 3: trapezoidal lattice h = 0.25 on [-24, 24]^2
    h, L = 0.25, 24.0
    m = int(L / h)
    tot = 0.0
    for i in range(-m, m + 1):
        for j in range(-m, m + 1):
            tot += S_of([i * h, j * h], 3)
    print("   N=3 lattice (h=%.2f): K = %.8f ; exact pi = %.8f" % (h, tot * h * h, math.pi))
    # N = 5: importance-sampled Monte Carlo in log coordinates, Laplace proposal of scale b
    b, M = 10.0, 40000
    acc, acc2 = 0.0, 0.0
    for _ in range(M):
        sv = [(-b if rng.random() < 0.5 else b) * math.log(1 - rng.random()) for _ in range(4)]
        p = 1.0
        for x in sv:
            p *= math.exp(-abs(x) / b) / (2 * b)
        val = S_of(sv, 5) / p
        acc += val
        acc2 += val * val
    mean = acc / M
    err = math.sqrt(max(acc2 / M - mean ** 2, 0) / M)
    print("   N=5 importance MC (%d points): K = %.3f +- %.3f ; exact pi^2 = %.3f" % (M, mean, err, math.pi ** 2))
    # (g) orientation of dr_2 ^ ... ^ dr_N against iota_E vol, and slice kernels with source
    print("\n(g) o_N = sgn of dr_2^...^dr_N relative to iota_E(dt_1...dt_N) (htk density convention, e0 = 0);")
    print("    slice t_red = s:  K_N(kappa; s) = o_N (-1)^n pi^{N/2}/Gamma(N/2) e^{s kappa^2} E[(sqrt(s) kappa + xi)^n]")
    for N in range(2, 8):
        n = N - 1
        osgn = set()
        for rep in range(3):
            t = tb.rand_t(N, rng)
            T = [Dual.var(t[e], e) for e in range(N)]
            tred = 1 / sum((1 / x for x in T), Dual(0))
            r = [tred / x for x in T]
            vol = {0: Fr(1)}
            for i in range(1, N):
                vol = wedge(vol, r[i].form())
            dens = htk_density(vol, t, N)
            osgn.add((dens > 0) - (dens < 0))
        o = osgn.pop() if len(osgn) == 1 else None
        poly = " + ".join("%s s^%d kappa^%d" % (Fr(math.comb(n, 2 * j)) * dfact(2 * j - 1) / Fr(2) ** j, (n - 2 * j) // 2 if (n - 2 * j) % 2 == 0 else 0, n - 2 * j)
                          for j in range(0, n // 2 + 1))
        print("   N=%d: o_N = %s;  E[(sqrt(s) k + xi)^n] = sum_j C(n,2j)(2j-1)!!/2^j (sqrt(s) k)^(n-2j)" % (N, o))
    print("\ntime %.1fs" % (time.time() - t0))


if __name__ == "__main__":
    main()
