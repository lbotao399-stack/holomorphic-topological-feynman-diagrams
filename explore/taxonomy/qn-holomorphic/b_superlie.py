"""Fork B: matrix Lie superalgebras gl(m|n), psl-free test algebras for 10d hCS (H=5, even supertrace pairing).
[E_ij, E_kl] = d_jk E_il - (-1)^{p(ij)p(kl)} d_li E_kj,  omega(E_ij, E_kl) = d_jk d_il (-1)^{|i|} (supertrace)."""
import itertools
import numpy as np
import superalg as sa


def glmn(m, n):
    N = m + n
    par_idx = [0] * m + [1] * n
    basis = [(i, j) for i in range(N) for j in range(N)]
    idx = {b: t for t, b in enumerate(basis)}
    D = len(basis)
    par = [(par_idx[i] + par_idx[j]) % 2 for i, j in basis]
    f = np.zeros((D, D, D), dtype=np.int64)
    for (i, j), A in idx.items():
        for (k, l), B in idx.items():
            if j == k:
                f[A, B, idx[(i, l)]] += 1
            if l == i:
                f[A, B, idx[(k, j)]] -= (-1) ** (par[A] * par[B])
    w = np.zeros((D, D), dtype=np.int64)
    for (i, j), A in idx.items():
        w[A, idx[(j, i)]] = (-1) ** par_idx[i]
    names = ['E%d%d' % b for b in basis]
    blocks = {'ev': [t for t in range(D) if par[t] == 0], 'od': [t for t in range(D) if par[t] == 1]}
    return sa.SuperLie('gl(%d|%d)' % (m, n), par, f, w, 0, names, blocks)


def str_ad_power(sl, k, rng, trials=3):
    """str over g_tot of (ad_X)^k for random even X (integer), exact."""
    out = []
    for _ in range(trials):
        ev = [t for t in range(sl.N) if sl.par[t] == 0]
        x = np.zeros(sl.N, dtype=object)
        for t in ev:
            x[t] = rng.randint(-3, 3)
        M = np.einsum('a,abc->cb', x, sl.f.astype(object))   # (ad X)[c,b]
        Mk = np.identity(sl.N, dtype=object)
        for _ in range(k):
            Mk = Mk.dot(M)
        s = sum(((-1) ** int(sl.par[t])) * Mk[t, t] for t in range(sl.N))
        out.append(s)
    return out


if __name__ == '__main__':
    import random
    rng = random.Random(3)
    for m, n in ((1, 1), (2, 1), (2, 2)):
        sl = glmn(m, n)
        sl.check()
        print(sl.name, 'dim', sl.N, 'check ok; str ad^6 on random even X:', str_ad_power(sl, 6, rng),
              'str ad^2:', str_ad_power(sl, 2, rng))
