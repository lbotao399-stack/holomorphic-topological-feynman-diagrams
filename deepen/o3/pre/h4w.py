import sys, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/deform/code')
import ce, la, spaces
w = int(sys.argv[1]); sl = sys.argv[2] == 'sl'
deg = (w, w)
t0 = time.time()
C2, k2 = spaces.basis(2, deg, sl)
C3, k3 = spaces.basis(3, deg, sl)
C4, k4 = spaces.basis(4, deg, sl)
print("weight", deg, "sl" if sl else "gl", "dim C2,C3,C4 =", len(C2), len(C3), len(C4), "t=%.1f" % (time.time()-t0), flush=True)
def rk(src, k, keys):
    if not src: return 0
    img = [spaces.const_part(ce.d_ce(spaces.to_np(c), k, sl)) for c in src]
    img = [x for x in img if x] 
    if not img: return 0
    _, M = la.matrix_from(img, keys)
    return la.certified_rank(M)['rank']
r2 = rk(C2, 2, k3); print("rank d2 =", r2, flush=True)
r3 = rk(C3, 3, k4); print("rank d3 =", r3, flush=True)
k5 = spaces.keys_of(5, deg, sl)
r4 = rk(C4, 4, k5); print("rank d4 =", r4, "t=%.1f" % (time.time()-t0), flush=True)
print("H3 =", len(C3) - r3 - r2, " H4 =", len(C4) - r4 - r3)
