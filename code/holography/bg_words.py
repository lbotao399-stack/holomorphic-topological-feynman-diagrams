"""Letter algebras (cobar word algebras A^!) for the tree-level twisted complexes.

A letter algebra has letters i = 0..n-1 with weight (tuple), parity (0/1), ghost number, and
Q_r on letters: Q[i] = list of (coeff, (a, b)) meaning Q_r x_i = sum coeff * x_a x_b.
Q_r acts on words as an odd derivation with Koszul signs.  The constant ghost is never a letter
(relative complex for the constant gauge algebra).

n4_letters(Lmax): superfield C(z, theta) = sum_A e_A X_A, e_A = z^m theta^S (S increasing),
    QC = C.C;  (-1)^{|S_C|} Q X_C = sum_{A.B=C} (-1)^{|X_A||e_B|} eps(S_A,S_B) X_A X_B,
    letter weight (theta1,theta2,theta3; n+, n-) = (1_S; m1, m2), parity 1+|S|, gh 1-|S|,
    level 2|S| + 3(m1+m2).
ctheta_letters(): A = C[theta]  (one even letter X, Q = 0).
cz_letters(Lmax):  A = C[z]      (odd letters c_n = coefficient of z^n, Q c_n = sum_{a+b=n} c_a c_b).
"""
from itertools import combinations


class LetterAlgebra:
    def __init__(self, names, weights, parities, ghs, Q):
        self.names = names
        self.w = weights
        self.par = parities
        self.gh = ghs
        self.Q = Q
        self.n = len(names)
        self.dimw = len(weights[0]) if weights else 0

    # ---- words
    def wword(self, word):
        s = [0] * self.dimw
        for i in word:
            for k, v in enumerate(self.w[i]):
                s[k] += v
        return tuple(s)

    def pword(self, word):
        return sum(self.par[i] for i in word) & 1

    def ghword(self, word):
        return sum(self.gh[i] for i in word)

    def Qword(self, word):
        """Q_r on a word: dict word -> coeff."""
        out = {}
        pre = 0
        for pos, i in enumerate(word):
            s = -1 if pre else 1
            for c, (a, b) in self.Q[i]:
                nw = word[:pos] + (a, b) + word[pos + 1:]
                out[nw] = out.get(nw, 0) + s * c
            pre ^= self.par[i]
        return {k: v for k, v in out.items() if v}

    # ---- cyclic words
    def canon_cyclic(self, word):
        """Tr(word) = sign * Tr(canonical rotation); sign = 0 if the trace vanishes identically."""
        n = len(word)
        if n == 0:
            return 1, word
        tot = self.pword(word)
        best = word
        best_sign = 1
        cur = word
        sign = 1
        seen_self_minus = False
        for k in range(1, n):
            x = cur[0]
            px = self.par[x]
            if px and (tot ^ px):
                sign = -sign
            cur = cur[1:] + (x,)
            if cur == word and sign == -1:
                seen_self_minus = True
            if cur < best:
                best = cur
                best_sign = sign
        if seen_self_minus:
            return 0, None
        return best_sign, best

    def words_upto(self, lam):
        """All nonempty words with weight <= lam componentwise (generator)."""
        lam = tuple(lam)
        letters = [i for i in range(self.n) if all(a <= b for a, b in zip(self.w[i], lam))]
        out = []

        def rec(word, rem):
            for i in letters:
                wi = self.w[i]
                if all(a <= b for a, b in zip(wi, rem)):
                    nw = word + (i,)
                    nrem = tuple(b - a for a, b in zip(wi, rem))
                    out.append(nw)
                    rec(nw, nrem)
        rec((), lam)
        return out


def _eps(SA, SB):
    """Sign of sorting the concatenation SA + SB (both increasing, disjoint)."""
    inv = 0
    for a in SA:
        for b in SB:
            if a > b:
                inv += 1
    return -1 if inv & 1 else 1


def n4_letters(Lmax):
    subsets = []
    for r in range(4):
        subsets += [tuple(c) for c in combinations((1, 2, 3), r)]
    keys = []
    for S in subsets:
        for m1 in range(Lmax // 3 + 1):
            for m2 in range(Lmax // 3 + 1 - m1):
                if not S and m1 == 0 and m2 == 0:
                    continue
                if 2 * len(S) + 3 * (m1 + m2) <= Lmax:
                    keys.append((m1, m2, S))
    keys.sort(key=lambda k: (2 * len(k[2]) + 3 * (k[0] + k[1]), k))
    index = {k: i for i, k in enumerate(keys)}
    names, weights, pars, ghs = [], [], [], []
    for (m1, m2, S) in keys:
        names.append(f"X[{m1},{m2};{''.join(map(str, S)) or '-'}]")
        weights.append(tuple(1 if I in S else 0 for I in (1, 2, 3)) + (m1, m2))
        pars.append((1 + len(S)) & 1)
        ghs.append(1 - len(S))
    Q = []
    for (m1, m2, S) in keys:
        terms = {}
        sC = -1 if len(S) & 1 else 1
        for r in range(len(S) + 1):
            for SA in combinations(S, r):
                SB = tuple(x for x in S if x not in SA)
                for a1 in range(m1 + 1):
                    for a2 in range(m2 + 1):
                        kA = (a1, a2, SA)
                        kB = (m1 - a1, m2 - a2, SB)
                        if kA not in index or kB not in index:
                            continue      # constant ghost (0,0,()) is not a letter
                        pXA = (1 + len(SA)) & 1
                        s = sC * (-1 if (pXA * len(SB)) & 1 else 1) * _eps(SA, SB)
                        key = (index[kA], index[kB])
                        terms[key] = terms.get(key, 0) + s
        Q.append([(c, ab) for ab, c in terms.items() if c])
    return LetterAlgebra(names, weights, pars, ghs, Q), keys


def ctheta_letters():
    return LetterAlgebra(["X"], [(1,)], [0], [0], [[]]), None


def cz_letters(Lmax):
    names = [f"c{n}" for n in range(1, Lmax + 1)]
    weights = [(n,) for n in range(1, Lmax + 1)]
    Q = []
    for n in range(1, Lmax + 1):
        Q.append([(1, (a - 1, n - a - 1)) for a in range(1, n)])
    return LetterAlgebra(names, weights, [1] * Lmax, [1] * Lmax, Q), None


def check_Q2(alg, maxlen_word=3):
    """Q_r^2 = 0 on every letter (as words)."""
    bad = []
    for i in range(alg.n):
        acc = {}
        for c, (a, b) in alg.Q[i]:
            for w, v in alg.Qword((a, b)).items():
                acc[w] = acc.get(w, 0) + c * v
        acc = {k: v for k, v in acc.items() if v}
        if acc:
            bad.append((alg.names[i], acc))
    return bad
