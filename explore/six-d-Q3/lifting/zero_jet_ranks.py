"""Independent finite Euler comparison for the zero-jet sl_N Koszul complex."""

from collections import Counter
from itertools import permutations
from math import factorial


def sign(p):
    return (-1) ** sum(i > j for k, i in enumerate(p) for j in p[k + 1:])


def adjoint_euler(n, cutoff):
    # Torus weights are integer n-tuples summing to zero.
    zero = (0,) * n
    roots = [tuple(int(k == i) - int(k == j) for k in range(n))
             for i in range(n) for j in range(n) if i != j]
    table = Counter({(0, 0, zero): 1})
    for weight in roots + [zero] * (n - 1):
        for dp, dq in ((1, 0), (0, 1), (1, 1)):
            updated = Counter()
            for key, value in table.items():
                p, q, w = key
                updated[key] += value
                power = 1
                while p + dp * power <= cutoff and q + dq * power <= cutoff:
                    ww = tuple(x + power * y for x, y in zip(w, weight))
                    odd = (dp, dq) == (1, 1)
                    updated[p + dp * power, q + dq * power, ww] += -value if odd else value
                    if odd:
                        break
                    power += 1
            table = {k: v for k, v in updated.items() if v}
    answer = Counter()
    for perm in permutations(range(n)):
        weight = tuple(i - perm[i] for i in range(n))
        for p in range(cutoff + 1):
            for q in range(cutoff + 1):
                answer[p, q] += sign(perm) * table.get((p, q, weight), 0)
    return answer


def cycles(p):
    sizes = []
    seen = set()
    for i in range(len(p)):
        if i in seen:
            continue
        size = 0
        while i not in seen:
            seen.add(i)
            i = p[i]
            size += 1
        sizes.append(size)
    return sizes


def multiply(a, b, limit):
    return [sum(a[i] * b[k - i] for i in range(k + 1)) for k in range(limit + 1)]


def cartan_euler(n, limit):
    answer = Counter()
    for perm in permutations(range(n)):
        denom = [1] + [0] * limit
        for k in cycles(perm):
            factor = [int(i % k == 0) for i in range(limit + 1)]
            denom = multiply(denom, factor, limit)
        # Remove the fixed center line of the permutation representation.
        sym = [v - (denom[i - 1] if i else 0) for i, v in enumerate(denom)]
        ext = [1] * (limit + 1)
        for k in cycles(perm):
            previous = ext
            ext = [v - (previous[i - k] if i >= k else 0) for i, v in enumerate(previous)]
        for p in range(limit + 1):
            for q in range(limit + 1):
                answer[p, q] += sum(ext[k] * sym[p - k] * sym[q - k]
                                    for k in range(min(p, q) + 1))
    for key, v in answer.items():
        assert v % factorial(n) == 0
        answer[key] = v // factorial(n)
    return answer


if __name__ == "__main__":
    for n, limit in ((2, 5), (3, 4), (4, 3)):
        a, b = adjoint_euler(n, limit), cartan_euler(n, limit)
        assert a == b, (n, a-b, b-a)
        print(f"sl_{n} adjoint and Cartan Euler coefficients agree in {(limit+1)**2} bidegrees")
