"""Murthy G_2^(3) to level 33 (pruned constant term) and the check Z_2 = Z_inf (1 + G1 + G2 + G3)."""
import sys, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
from holography import giant as g
L = int(sys.argv[1]) if len(sys.argv) > 1 else 33
N = 2
t = time.time()
G1 = g.G1(N, L); G2 = g.G2(N, L)
G3 = g.G3(N, L)
print("G3 time %.0fs" % (time.time() - t))
Zi = g.Z_inf(L); Z2 = g.Z_U(2, L)
S = [(1 if j == 0 else 0) + G1[j] + G2[j] + G3[j] for j in range(L + 1)]
lhs = g.ser_mul1(Zi, S, L)
print("G1", G1); print("G2", G2); print("G3", G3)
print("Z_U(2)            ", Z2)
print("Zinf(1+G1+G2+G3)  ", lhs)
print("diff", [a - b for a, b in zip(Z2, lhs)])
