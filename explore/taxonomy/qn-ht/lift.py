"""d_2 = [Q_2^hat] on H(Q_0), sector by sector (exact Fractions).
Sector: fixed torus weight vector w (letter weight = spin + alpha) and conserved charges (charge 0 for
the gauge Cartan); complex graded by ghost number (Q raises it by gstep) or by parity (Z/2).
Usage: python3 lift.py <spec> <wmax> [charge-filter]   Log: lift_<spec>_<wmax>.log
wmax: sectors with sum(w) <= wmax (weights in units of 1/2 for spin-1/2 letters: integers after
scaling by 2), w non-increasing (Weyl group of U(H)), only gauge-neutral charge vectors."""
import sys, time, itertools
from fractions import Fraction as Fr
sys.path.insert(0, ".")
import hteng
from specs import SPECS


class Ech:
    def __init__(self):
        self.rows = {}

    def reduce(self, v):
        v = dict(v)
        while v:
            p = min(v)
            if p in self.rows:
                c = v[p]
                for k, x in self.rows[p].items():
                    nv = v.get(k, 0) - c * x
                    if nv == 0:
                        v.pop(k, None)
                    else:
                        v[k] = nv
            else:
                return v
        return v

    def add(self, v):
        v = self.reduce(v)
        if not v:
            return False
        p = min(v)
        c = v[p]
        self.rows[p] = {k: x / c for k, x in v.items()}
        return True

    def rank(self):
        return len(self.rows)


def kernel_basis(images):
    piv, ker = {}, []
    for j, img in enumerate(images):
        v = dict(img)
        comb = {j: Fr(1)}
        while v:
            p = min(v)
            if p in piv:
                pi, pc = piv[p]
                c = v[p]
                for k, x in pi.items():
                    nv = v.get(k, 0) - c * x
                    if nv == 0:
                        v.pop(k, None)
                    else:
                        v[k] = nv
                for k, x in pc.items():
                    nv = comb.get(k, 0) - c * x
                    if nv == 0:
                        comb.pop(k, None)
                    else:
                        comb[k] = nv
            else:
                break
        if not v:
            ker.append(comb)
        else:
            p = min(v)
            c = v[p]
            piv[p] = ({k: x / c for k, x in v.items()}, {k: x / c for k, x in comb.items()})
    return ker


class Lifter:
    def __init__(self, specname, seed=1):
        self.sp = sp = SPECS[specname]
        self.th = th = sp['build']()
        self.H = th.H
        self.E = hteng.Engine(th, seed=seed)
        self.nq2nz = 0  # cocycle basis vectors z with Q2^hat z != 0 (as cochains), target H != 0
        # integer weights: scale by 2 if any half-integer spin
        self.base = th.letter_pool()
        from math import lcm
        self.scale = 1
        for L in self.base:
            for x in sp['spin'](L[0], L[1]):
                self.scale = lcm(self.scale, Fr(x).denominator)

    def weight(self, L):
        s = self.sp['spin'](L[0], L[1])
        return tuple(int(self.scale * (Fr(s[i]) + L[2][i])) for i in range(self.H))

    def grade(self, m):
        g = self.sp['ghost']
        if all(g(L[0], L[1]) is not None for L in m) and self.sp.get('gstep'):
            return sum(g(L[0], L[1]) for L in m)
        return sum(self.th.parity(L) for L in m) % 2

    def zgraded(self):
        return all(self.sp['ghost'](L[0], L[1]) is not None for L in self.base)

    def charge(self, m):
        ch = self.sp['charge']
        tot = None
        for L in m:
            c = ch(L[0], L[1])
            tot = c if tot is None else tuple(a + b for a, b in zip(tot, c))
        return tot

    def letters(self, wmax):
        out = []
        for (n, i, _) in self.base:
            for a in itertools.product(range(wmax + 1), repeat=self.H):
                L = (n, i, a)
                w = self.weight(L)
                if sum(w) <= wmax:
                    out.append(L)
        return sorted(out)

    def sector(self, w, letters):
        """monomials of total weight exactly w (tuple); grouped {charge: {grade: [monos]}}"""
        th = self.th
        Ls = [(L, self.weight(L)) for L in letters if all(a <= b for a, b in zip(self.weight(L), w))]
        res = {}
        H = self.H

        def rec(start, cur, rem):
            if not any(rem):
                s, m = th.canon(cur)
                if s:
                    ch = self.charge(m) if m else None
                    res.setdefault(ch, {}).setdefault(self.grade(m), set()).add(m)
            for idx in range(start, len(Ls)):
                L, wl = Ls[idx]
                if any(a > b for a, b in zip(wl, rem)):
                    continue
                if not any(wl):
                    # weight-0 letters (c^a): allowed, odd, no repetition
                    pass
                if th.parity(L) and cur and cur[-1] == L:
                    continue
                rec(idx, cur + [L], tuple(b - a for a, b in zip(wl, rem)))
        # weight-0 letters only finitely many (odd) - guard: they are odd so no repetition
        rec(0, [], w)
        return {ch: {g: sorted(v) for g, v in d.items()} for ch, d in res.items()}

    def run_sector(self, w, letters, out, show=2, chfilter=None):
        E = self.E
        secs = self.sector(w, letters)
        step = self.sp.get('gstep', 1) if self.zgraded() else 1
        results = []
        for ch in sorted(secs, key=lambda x: (x is None, x)):
            if ch is None:
                continue
            if chfilter is not None and not chfilter(ch):
                continue
            bas = secs[ch]
            zg = self.zgraded()
            def nxt(g):
                return g + step if zg else (g + 1) % 2
            idx = {g: {m: i for i, m in enumerate(b)} for g, b in bas.items()}
            Z, Bim = {}, {}
            for g, b in bas.items():
                imgs = []
                for m in b:
                    P = E.Q0({m: Fr(1)})
                    tgt = idx.get(nxt(g), {})
                    v = {}
                    for mm, c in P.items():
                        v[tgt[mm]] = c
                    imgs.append(v)
                Z[g] = kernel_basis(imgs)
                ech = Ech()
                for v in imgs:
                    if v:
                        ech.add(v)
                Bim[nxt(g)] = ech
            Hd = {g: len(Z[g]) - (Bim[g].rank() if g in Bim else 0) for g in bas}
            Hn = {g: h for g, h in sorted(Hd.items()) if h}
            if not Hn:
                continue
            ranks, lifted = {}, []
            nq0 = self.nq2nz
            for g, h in Hn.items():
                g2 = nxt(g)
                if not Hd.get(g2):
                    continue
                base = Bim.get(g2, Ech())
                ech = Ech()
                ech.rows = dict(base.rows)
                r0 = ech.rank()
                for zv in Z[g]:
                    P = {bas[g][i]: c for i, c in zv.items()}
                    Q2P = E.Q2(P)
                    if not Q2P:
                        continue
                    self.nq2nz += 1
                    v = {idx[g2][mm]: c for mm, c in Q2P.items()}
                    if ech.add(v):
                        lifted.append((g, P, Q2P))
                if ech.rank() > r0:
                    ranks[g] = ech.rank() - r0
            line = "w=%s charges=%s: dim C = %s ; dim H(Q0) = %s ; Q2^hat z != 0 on %d cocycle basis vectors%s" % (
                tuple(Fr(x, self.scale) for x in w) if self.scale != 1 else w, ch,
                {g: len(b) for g, b in sorted(bas.items())}, Hn, self.nq2nz - nq0,
                (" ; rank d2 = %s" % ranks) if ranks else "")
            out(line)
            for (g, P, Q2P) in lifted[:show]:
                out("    lifted (grade %s): O = %s" % (g, hteng.pstr(P, 12)))
                out("       Q2^hat O = %s" % hteng.pstr(Q2P, 12))
            results.append((w, ch, Hn, ranks))
        return results


def weights_upto(H, wmax):
    out = []
    for w in itertools.product(range(wmax + 1), repeat=H):
        if sum(w) <= wmax and all(w[i] >= w[i + 1] for i in range(H - 1)) and sum(w) > 0:
            out.append(w)
    return sorted(out, key=lambda w: (sum(w), w))


if __name__ == "__main__":
    name, wmax = sys.argv[1], int(sys.argv[2])
    log = open("lift_%s_%d.log" % (name, wmax), "w")
    def out(*a):
        s = " ".join(str(x) for x in a)
        print(s); log.write(s + "\n"); log.flush()
    t0 = time.time()
    Lf = Lifter(name)
    out("theory:", Lf.th.name, "(H,T)=(%d,1); grading:" % Lf.H, "ghost" if Lf.zgraded() else "parity",
        "; weights in units of 1/%d" % Lf.scale, "; sectors with sum(w) <= %d, w non-increasing, gauge charge 0" % wmax)
    letters = Lf.letters(wmax)
    out("letters:", len(letters))
    for w in weights_upto(Lf.H, wmax):
        Lf.run_sector(w, letters, out, chfilter=lambda ch: ch[0] == 0)
    out("time %.0fs" % (time.time() - t0))
