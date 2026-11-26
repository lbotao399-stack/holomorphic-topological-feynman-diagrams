"""Per-sector index residual of the k=1 truncation for U(2):
   res(lam) = [Z_U(2)]_lam - [Z_inf (1 + G_2^(1))]_lam  (= [Z_inf G_2^(2)]_lam below level 30).
usage: python3 a1_residual.py Jmin Jmax"""
import sys, time
from multiprocessing import Pool
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
from holography import giant as gi

def reps(j):
    out = []
    for n in range(j // 3 + 1):
        r = j - 3 * n
        if r % 2:
            continue
        Th = r // 2
        for npl in range(n, (n + 1) // 2 - 1, -1):
            nmi = n - npl
            if npl < nmi:
                continue
            for t1 in range(Th, -1, -1):
                for t2 in range(min(t1, Th - t1), -1, -1):
                    t3 = Th - t1 - t2
                    if t3 > t2 or t3 < 0:
                        continue
                    out.append((t1, t2, t3, npl, nmi))
    return out

def work(lam):
    t0 = time.time()
    zu = gi.Z_U2_box(lam).get(lam, 0)
    zi = gi.Zinf_box(lam)
    g1 = gi.G1_box(2, lam)
    one = dict(g1); z0 = (0,)*5
    one[z0] = one.get(z0, 0) + 1
    tot = 0
    for e1, c1 in zi.items():
        e2 = tuple(a - b for a, b in zip(lam, e1))
        if min(e2) < 0: continue
        tot += c1 * one.get(e2, 0)
    return lam, zu, tot, zu - tot, time.time() - t0

if __name__ == '__main__':
    J0, J1 = int(sys.argv[1]), int(sys.argv[2])
    sectors = [lam for j in range(J0, J1 + 1) for lam in reps(j)]
    print('# level | sector | Z_U2 | Z_inf(1+G1) | residual | time', flush=True)
    with Pool(12) as p:
        for lam, zu, tot, res, dt in p.imap(work, sectors):
            lev = 2*sum(lam[:3]) + 3*sum(lam[3:])
            print(f"{lev:3d} | {lam} | {zu} | {tot} | {res} | {dt:.1f}s", flush=True)
