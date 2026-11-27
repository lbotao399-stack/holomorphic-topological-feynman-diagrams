"""Polynomials over F_p (coefficient lists, low degree first): interpolation, evaluation, number of roots."""
# ---- polynomials mod p (lists of coefficients, low degree first)
def pmod(a, p):
    a = [x % p for x in a]
    while a and a[-1] == 0:
        a.pop()
    return a


def pdivmod(a, b, p):
    a = pmod(a, p); b = pmod(b, p)
    q = [0] * max(len(a) - len(b) + 1, 1)
    inv = pow(b[-1], p - 2, p)
    while len(a) >= len(b) and a:
        c = a[-1] * inv % p
        k = len(a) - len(b)
        q[k] = c
        for i, x in enumerate(b):
            a[i + k] = (a[i + k] - c * x) % p
        a = pmod(a, p)
    return pmod(q, p), a


def pmul(a, b, p):
    out = [0] * (len(a) + len(b) - 1) if a and b else []
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] = (out[i + j] + x * y) % p
    return pmod(out, p)


def pgcd(a, b, p):
    a, b = pmod(a, p), pmod(b, p)
    while b:
        a, b = b, pdivmod(a, b, p)[1]
    return a


def n_roots_mod(f, p):
    """number of distinct roots of f in F_p: deg gcd(f, x^p - x)."""
    f = pmod(f, p)
    if len(f) <= 1:
        return 0 if f else None
    # x^p mod f by square and multiply
    res, base, e = [1], [0, 1], p
    while e:
        if e & 1:
            res = pdivmod(pmul(res, base, p), f, p)[1]
        base = pdivmod(pmul(base, base, p), f, p)[1]
        e >>= 1
    xpx = pmod([(res[i] if i < len(res) else 0) - (1 if i == 1 else 0) for i in range(max(len(res), 2))], p)
    g = pgcd(f, xpx, p)
    return len(g) - 1


def interp(xs, ys, p):
    """Lagrange interpolation mod p -> coefficient list."""
    n = len(xs)
    coef = [0] * n
    for i in range(n):
        num, den = [1], 1
        for j in range(n):
            if j != i:
                num = pmul(num, [(-xs[j]) % p, 1], p)
                den = den * (xs[i] - xs[j]) % p
        c = ys[i] * pow(den, p - 2, p) % p
        for k, x in enumerate(num):
            coef[k] = (coef[k] + c * x) % p
    return pmod(coef, p)


def peval(c, x, p):
    return sum(ci * pow(x, i, p) for i, ci in enumerate(c)) % p


