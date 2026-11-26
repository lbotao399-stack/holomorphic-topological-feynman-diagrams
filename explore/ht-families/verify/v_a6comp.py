"""a6, second part: number of nonzero coordinate components of rho_4 in the engine's small-plan ratio
coordinates (own Pf, own cycle basis, own exact Jacobian of the engine t-map), fan G3 and G4; also large plan."""
import sys, random
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/ht-families/verify")
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
import vcore as vc
import htk
from htk import spred
from v_b1 import GR
rng = random.Random(4)
for k in (3, 4):
    for prefer in ("small", "large"):
        G = htk.Graph(6, GR[k])
        red = spred.sp_reduce(G, prefer=prefer)
        net = vc.Net(6, GR[k])
        C = net.cycle_basis()
        sizes = []
        for _ in range(3):
            pt = red.complete_point({c: v for c, v in red.random_point(rng).items() if c in red.coords})
            t = red.tmap(pt)
            dts = []
            for e in range(9):
                vec = []
                for c in red.coords:
                    blk = next(b for b in red.blocks if c in b.bary)
                    b0 = blk.bary[0]
                    vec.append(t[e] * (Fr(red.texp[e].get(c, 0)) / pt[c] - Fr(red.texp[e].get(b0, 0)) / pt[b0]))
                dts.append(vc.oneform(vec))
            pf, U, _ = vc.rho_forms(C, t, dts)
            sizes.append(len(pf))
        print("G%d plan %-5s coordinates %s: #nonzero components of rho_4 at 3 rational points: %s" % (k, prefer, red.coords, sizes))
