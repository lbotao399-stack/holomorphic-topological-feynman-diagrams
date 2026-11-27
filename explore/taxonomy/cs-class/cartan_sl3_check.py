"""sl_3 consistency of small_N (o' = -(1/256) Q2323 + d eta) with the Cartan detector:
o'|_C = 81 det(h (x) p) (cartan_detector.log), hence Q2323|_C must be -20736 det(h (x) p)."""
import numpy as np
from fractions import Fraction as Fr
from common import Log, PRIMES, det_int, ratrec, crt_pair
import lie
import cartan_lib as cl

log = Log(__file__)
L = lie.build('A2')
W = [[(1 if k == i else 0) - (1 if k == i + 1 else 0) for i in range(2)] for k in range(3)]
P23 = lambda x: sum(sum(w[i] * x[i] for i in range(2)) ** 2 for w in W) * sum(sum(w[i] * x[i] for i in range(2)) ** 3 for w in W)
res = {}
for p in PRIMES[:2]:
    T = cl.sym_tensor_from_poly(P23, 2, p)
    det_ = cl.Detector(L, p)
    rng = np.random.default_rng(5)
    out = []
    for t in range(6):
        hc, mom = det_.sample(rng)
        kd = det_int([[int(hc[v, i]) * int(mom[v, k]) for i in range(2) for k in range(4)] for v in range(8)]) % p
        if kd == 0:
            out.append(None); continue
        q = cl.half_bracket_tensor(L, T, hc, mom, p)
        o = det_.o_prime(hc, mom)
        out.append((q * pow(kd, p - 2, p) % p, o * pow(kd, p - 2, p) % p))
    res[p] = out
for a, b in zip(res[PRIMES[0]], res[PRIMES[1]]):
    if a is None or b is None:
        log('degenerate sample (det = 0) skipped'); continue
    (q0, o0), (q1, o1) = a, b
    rq = ratrec(*crt_pair(q0, PRIMES[0], q1, PRIMES[1]))
    ro = ratrec(*crt_pair(o0, PRIMES[0], o1, PRIMES[1]))
    log('Q2323|_C / det = %s,  o\'|_C / det = %s,  o\'/Q2323 on Cartan = %s' % (rq, ro, ro / rq))
log.close()
