"""S2 eq-IA, eq-IB (S2 labelling) evaluated on sum lambda = 0, against Delta' = det(l1, l2, l4) (A),
Delta'' = det(lu, lw, la) (B); and K of S2 graph A (own GM integration) at the same sources."""
import sys, random
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/anomaly-structure/verify')
import vk, sp_coords as spc
rng = random.Random(5)
det3 = lambda a, b, c: vk.det([[a[i], b[i], c[i]] for i in range(3)])
for _ in range(3):
    l1, l2, l4 = [[rng.randint(-6, 6) for _ in range(3)] for _ in range(3)]
    l5 = [-(l1[i] + l2[i] + l4[i]) for i in range(3)]
    L = {1: l1, 2: l2, 4: l4, 5: l5}
    b = lambda i, j, k: det3(L[i], L[j], L[k])
    IA = Fr(1, 3) * (b(1, 2, 4) / Fr(12) + b(1, 2, 5) / Fr(24)) * (b(1, 4, 5) / Fr(6) + b(2, 4, 5) / Fr(3))
    Dp = b(1, 2, 4)
    K = spc.integrate_paths(6, [(5, 0), (0, 1), (1, 2), (5, 2), (2, 3), (3, 4), (4, 5)], [[0, 1, 2], [3], [6, 5, 4]],
                            {0: l1, 1: l2, 2: [0, 0, 0], 3: l4, 4: l5}, 3, 1)
    print("A: Delta' = %s, I_A(eq-IA) = %s, -Delta'^2/432 = %s, K(S2 labelling) = %s" % (Dp, IA, -Fr(Dp ** 2, 432), K))
    lu, lw, la = [[rng.randint(-6, 6) for _ in range(3)] for _ in range(3)]
    lb = [-(lu[i] + lw[i] + la[i]) for i in range(3)]
    L = {'u': lu, 'w': lw, 'a': la, 'b': lb}
    b = lambda i, j, k: det3(L[i], L[j], L[k])
    IB = -Fr(1, 6) * (b('u', 'w', 'a') / Fr(12) + b('u', 'w', 'b') / Fr(24)) * (b('u', 'a', 'b') / Fr(2) + b('w', 'a', 'b') / Fr(2))
    print("B: I_B(eq-IB) on sum = 0: %s" % IB)
