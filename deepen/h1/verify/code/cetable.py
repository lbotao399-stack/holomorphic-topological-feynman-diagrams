"""H^n(C_g^0)_(a,b) from the verifier's jet code (cost (0,0,0,a,b) slice of V), mod P1, P2 (and Q if 'exact')."""
import sys, time
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from vjet import Jets, P1, P2
from vcohom import cohom_V
alg = sys.argv[1]; charge0 = 'charge0' in sys.argv; exact = 'exact' in sys.argv
abl = [tuple(int(x) for x in s.split(',')) for s in sys.argv[2].split(';')]
J = Jets(alg)
for (a, b) in abl:
    t0 = time.time()
    out, S = cohom_V(J, (0, 0, 0, a, b), charge0, exact=exact)
    print(f"{alg} charge0={charge0} (a,b)=({a},{b}) " + " ".join(f"n={g}:dimC={o['dimM']},dimC0={o['dimV']},H={[o[P1], o[P2]] + ([o[None]] if exact else [])}" for g, o in sorted(out.items())), f"[{time.time()-t0:.1f}s]", flush=True)
