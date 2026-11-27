"""Fork B: runtime patch of htk.spred._poly_det: for n >= 8 use the column-subset dynamic programme (as for n < 8)
instead of Bareiss with exact sympy division (python-flint absent; Bareiss on the 8x8 local matrix of the 3-chain
move of theta(3,3,3) with 24 symbolic output momenta did not finish in 5 min).  Import before building densities."""
import htk.spred as _S
_orig = _S._poly_det


def poly_det_dp(rows, ring):
    n = len(rows)
    if n < 8:
        return _orig(rows, ring)
    order = sorted(range(n), key=lambda i: sum(1 for x in rows[i] if x != 0))
    perm = list(order)
    sgn = 1
    for i in range(n):
        while perm[i] != i:
            j = perm[i]
            perm[i], perm[j] = perm[j], perm[i]
            sgn = -sgn
    dp = {0: ring.one}
    for i in order:
        row = rows[i]
        nz = [c for c in range(n) if row[c] != 0]
        new = {}
        for mask, val in dp.items():
            for c in nz:
                bit = 1 << c
                if mask & bit:
                    continue
                term = row[c] * val
                if bin(mask >> (c + 1)).count("1") & 1:
                    term = -term
                m2 = mask | bit
                new[m2] = new[m2] + term if m2 in new else term
        dp = {m: v for m, v in new.items() if v != 0}
        if not dp:
            return ring.zero
    res = dp.get((1 << n) - 1, ring.zero)
    return res if sgn > 0 else -res


_S._poly_det = poly_det_dp
