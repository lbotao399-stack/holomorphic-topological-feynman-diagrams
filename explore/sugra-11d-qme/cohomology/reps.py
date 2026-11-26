"""Exact cocycle representatives for H(C_quot) and H(V) on a slot (w, t, N).
H(C): ker(Q: C^N -> C^{N+1}) reduced modulo Q(C^{N-1}).
H(V): {z in C^N : Qz in I^{N+1}} reduced modulo Q(C^{N-1}) + I^N,  I = sum_k D_k C (slot t + e_k).
Kernel vectors by integer elimination with tag columns; the returned representatives are fully reduced against
the echelon of the subspace that is divided out (canonical for the fixed monomial order)."""
from math import gcd
from jets import Echelon, vec


def _kernel(vectors, pre=None):
    """vectors: list of (real dict, tag dict).  pre: Echelon of real vectors already present (no tags).
    returns integer tag-combinations z with sum z_j vectors_j in span(pre)."""
    piv = {}  # real lead -> (real, tag)
    prepiv = pre.piv if pre is not None else {}
    ker = []
    for real, tag in vectors:
        real = dict(real)
        tag = dict(tag)
        while real:
            c = max(real)
            if c in prepiv:
                p = prepiv[c]
                a, b = p[c], real[c]
                g = gcd(a, b); a //= g; b //= g
                real = {k: a * x for k, x in real.items()}
                tag = {k: a * x for k, x in tag.items()}
                for k, x in p.items():
                    y = real.get(k, 0) - b * x
                    if y:
                        real[k] = y
                    else:
                        real.pop(k, None)
                continue
            if c in piv:
                pr, pt = piv[c]
                a, b = pr[c], real[c]
                g = gcd(a, b); a //= g; b //= g
                real = {k: a * x for k, x in real.items()}
                tag = {k: a * x for k, x in tag.items()}
                for k, x in pr.items():
                    y = real.get(k, 0) - b * x
                    if y:
                        real[k] = y
                    else:
                        real.pop(k, None)
                for k, x in pt.items():
                    y = tag.get(k, 0) - b * x
                    if y:
                        tag[k] = y
                    else:
                        tag.pop(k, None)
                continue
            piv[c] = (real, tag)
            break
        if not real:
            g = 0
            for x in tag.values():
                g = gcd(g, x)
            ker.append({k: x // g for k, x in tag.items()})
    return ker


def full_reduce(E, v):
    """reduce v completely (all entries, not only the lead) against the echelon E; integer-scaled result."""
    v = {k: x for k, x in v.items() if x}
    out = {}
    while v:
        c = max(v)
        p = E.piv.get(c)
        if p is None:
            out[c] = v.pop(c)
            continue
        a, b = p[c], v[c]
        g = gcd(a, b); a //= g; b //= g
        v = {k: a * x for k, x in v.items()}
        out = {k: a * x for k, x in out.items()}
        for k, x in p.items():
            y = v.get(k, 0) - b * x
            if y:
                v[k] = y
            else:
                v.pop(k, None)
    g = 0
    for x in out.values():
        g = gcd(g, x)
    if g > 1:
        out = {k: x // g for k, x in out.items()}
    return out


def reps_C(Co, t, N):
    A = Co.A
    mid = Co.basis(t, N)
    tgt = Co.basis(t, N + 1)
    src = Co.basis(t, N - 1)
    tidx = {m: j for j, m in enumerate(tgt)}
    midx = {m: j for j, m in enumerate(mid)}
    vecs = [(vec(A.Qmono(m), tidx), {j: 1}) for j, m in enumerate(mid)]
    ker = _kernel(vecs)
    E = Echelon()
    for m in src:
        E.insert(vec(A.Qmono(m), midx))
    reps = []
    for z in ker:
        if E.insert(dict(z)):
            reps.append(z)
    # canonical: reduce each representative modulo the image and the other representatives' leads
    Eim = Echelon()
    for m in src:
        Eim.insert(vec(A.Qmono(m), midx))
    out = []
    for z in reps:
        r = full_reduce(Eim, z)
        out.append({mid[j]: c for j, c in r.items()})
    return out


def reps_V(Co, t, N):
    A = Co.A
    n = A.n
    mid = Co.basis(t, N)
    tgt = Co.basis(t, N + 1)
    src = Co.basis(t, N - 1)
    tidx = {m: j for j, m in enumerate(tgt)}
    midx = {m: j for j, m in enumerate(mid)}
    Itgt = Echelon()
    for k in range(n):
        for m in Co.basis(Co.shiftw(t, k), N + 1):
            Itgt.insert(vec(A.Dq({m: 1}, k), tidx))
    vecs = [(vec(A.Qmono(m), tidx), {j: 1}) for j, m in enumerate(mid)]
    ker = _kernel(vecs, pre=Itgt)
    E = Echelon()
    for k in range(n):
        for m in Co.basis(Co.shiftw(t, k), N):
            E.insert(vec(A.Dq({m: 1}, k), midx))
    for m in src:
        E.insert(vec(A.Qmono(m), midx))
    reps = []
    Ecopy = Echelon()
    Ecopy.piv = {k: dict(v) for k, v in E.piv.items()}
    for z in ker:
        if E.insert(dict(z)):
            reps.append(z)
    out = []
    for z in reps:
        r = full_reduce(Ecopy, z)
        out.append({mid[j]: c for j, c in r.items()})
    return out
