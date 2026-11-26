"""PSLQ tests on a high-precision value of c (mpmath)."""
import sys
import mpmath as mp

def tests(cval, dps, tol_digits, maxcoeff=10**7):
    mp.mp.dps = dps
    z3, pi2, l2 = mp.zeta(3), mp.pi ** 2, mp.log(2)
    tol = mp.mpf(10) ** (-tol_digits)
    out = {}
    for name, basis in [("1", [cval, 1]),
                        ("1,z3", [cval, 1, z3]),
                        ("1,z3,pi2,log2", [cval, 1, z3, pi2, l2])]:
        rel = mp.pslq(basis, tol=tol, maxcoeff=maxcoeff, maxsteps=10**6)
        out[name] = rel
    return out

if __name__ == "__main__":
    mp.mp.dps = 80
    cval = mp.mpf(sys.argv[1])
    print(tests(cval, 80, int(sys.argv[2])))
