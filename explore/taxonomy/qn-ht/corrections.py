"""Lowest H(Q_0) representatives z with Q2^hat z != 0 (as a cochain), and the two-loop correction:
Q2^hat z = Q_0 y  =>  Q(z - N_2 y) = O(hbar^4)  (no lift).  Scans sectors in increasing weight.
Usage: python3 corrections.py <spec> <wmax> [nshow]   Log: corrections_<spec>.log"""
import sys, time
from fractions import Fraction as Fr
sys.path.insert(0, ".")
import hteng
from lift import Lifter, Ech, kernel_basis, weights_upto


def solve(images, v):
    """find x with sum_j x_j images[j] = v (dicts over target index); None if inconsistent"""
    piv = {}
    for j, img in enumerate(images):
        r = dict(img)
        comb = {j: Fr(1)}
        while r:
            p = min(r)
            if p in piv:
                pr, pc = piv[p]
                c = r[p]
                for k, x in pr.items():
                    nv = r.get(k, 0) - c * x
                    if nv == 0:
                        r.pop(k, None)
                    else:
                        r[k] = nv
                for k, x in pc.items():
                    nv = comb.get(k, 0) - c * x
                    if nv == 0:
                        comb.pop(k, None)
                    else:
                        comb[k] = nv
            else:
                c = r[p]
                piv[p] = ({k: x / c for k, x in r.items()}, {k: x / c for k, x in comb.items()})
                break
    r = dict(v)
    sol = {}
    while r:
        p = min(r)
        if p not in piv:
            return None
        pr, pc = piv[p]
        c = r[p]
        for k, x in pr.items():
            nv = r.get(k, 0) - c * x
            if nv == 0:
                r.pop(k, None)
            else:
                r[k] = nv
        for k, x in pc.items():
            sol[k] = sol.get(k, 0) + c * x
    return {k: x for k, x in sol.items() if x}


name, wmax = sys.argv[1], int(sys.argv[2])
nshow = int(sys.argv[3]) if len(sys.argv) > 3 else 1
log = open("corrections_%s.log" % name, "w")
def out(*a):
    s = " ".join(str(x) for x in a)
    print(s); log.write(s + "\n"); log.flush()
t0 = time.time()
Lf = Lifter(name)
E, th = Lf.E, Lf.th
out("theory:", th.name, "; weights in units of 1/%d; first %d sectors with a representative z, Q2^hat z != 0" % (Lf.scale, nshow))
letters = Lf.letters(wmax)
found = 0
for w in weights_upto(Lf.H, wmax):
    secs = Lf.sector(w, letters)
    zg = Lf.zgraded()
    step = Lf.sp.get('gstep', 1) if zg else 1
    nxt = (lambda g: g + step) if zg else (lambda g: (g + 1) % 2)
    prv = (lambda g: g - step) if zg else (lambda g: (g + 1) % 2)
    for ch in sorted(k for k in secs if k is not None and k[0] == 0):
        bas = secs[ch]
        idx = {g: {m: i for i, m in enumerate(b)} for g, b in bas.items()}
        def q0img(g):
            tgt = idx.get(nxt(g), {})
            return [{tgt[mm]: c for mm, c in E.Q0({m: Fr(1)}).items()} for m in bas[g]]
        for g in sorted(bas):
            imgs = q0img(g)
            Z = kernel_basis(imgs)
            if not Z:
                continue
            ech = Ech()
            if prv(g) in bas:
                for v in q0img(prv(g)):
                    if v:
                        ech.add(v)
            for zv in Z:
                if not ech.add(dict(zv)):
                    continue  # exact (or dependent)
                P = {bas[g][i]: c for i, c in zv.items()}
                Q2P = E.Q2(P)
                if not Q2P:
                    continue
                v = {idx[nxt(g)][mm]: c for mm, c in Q2P.items()}
                # solve Q0 y = Q2P with y in C^g
                y = solve(imgs, v)
                wl = tuple(Fr(x, Lf.scale) for x in w)
                out("w=%s charges=%s grade %s: z = %s" % (wl, ch, g, hteng.pstr(P, 14)))
                out("    Q2^hat z = %s" % hteng.pstr(Q2P, 14))
                if y is None:
                    out("    NOT Q0-exact (lift)")
                else:
                    Y = {bas[g][i]: c for i, c in y.items()}
                    chk = hteng.psum(E.Q0(Y), {m: -c for m, c in Q2P.items()})
                    out("    y = %s ;  Q0 y - Q2^hat z = %s" % (hteng.pstr(Y, 14), hteng.pstr(chk)))
                found += 1
                break
            if found >= nshow:
                break
        if found >= nshow:
            break
    if found >= nshow:
        break
out("[%.0fs]" % (time.time() - t0))
