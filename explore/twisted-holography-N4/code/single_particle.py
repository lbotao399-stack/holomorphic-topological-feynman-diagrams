"""Twisted IIB on X = C^5 \\ C^2 (Costello-Li 1606.00365 Sec. 9) vs planar single-trace classes of the
holomorphic twist of N=4 U(N) SYM.  Exact ranks over Q (fractions).

Bulk (t^0, Dolbeault cohomology of X):  A0 = C[z1,z2,w1,w2,w3] (x) Lambda(xi),  A2 = C[z] (x) L (x) Lambda(xi),
  L = span{ w^(-a-1) : a in N^3 } = H^2(C^3\\0, O),  xi = (xi_z1, xi_z2, xi_w1, xi_w2, xi_w3) odd (polyvectors d/dx).
  Divergence  del = sum_x d/dx d/dxi_x  (Omega = dz1 dz2 dw1 dw2 dw3).
  F = f xi_z1 xi_z2, f = (w1 w2 w3)^(-1).  For del(alpha)=0: {F, alpha} = del(F alpha);  F . w^b = 0 unless b = 0.
  coker(W) = dim ker(del | A2_W) - rank del( F . ker(del | A0_{W - w(F)}) ),   w(F) = (-1,-1,-1,-1,-1).
Boundary: de Rham complex of C^{2|3}: z (even), theta (odd), dz (odd), dtheta (even).
  Prediction:  coker(W) = dim B^{>=2}(beta(W)) + [beta(W) = 0],   beta(W) = (W_z + 1, -W_w - 1)  (n_a, theta_I).
"""
from fractions import Fraction
from itertools import product, combinations
import sys

def rank(rows):
    rows = [dict(r) for r in rows if r]
    piv = {}
    r = 0
    for row in rows:
        row = {k: Fraction(v) for k, v in row.items() if v}
        for k in sorted(row):
            pass
        # reduce
        changed = True
        while row:
            k = min(row)
            if k in piv:
                prow = piv[k]
                c = row[k] / prow[k]
                for kk, vv in prow.items():
                    nv = row.get(kk, 0) - c * vv
                    if nv: row[kk] = nv
                    else: row.pop(kk, None)
            else:
                piv[k] = row
                r += 1
                break
    return r

def kernel_basis(cols, ncols_keys):
    """cols: list of (vector dict) images of basis elements; return basis of kernel as list of dicts over basis idx."""
    # Gaussian elimination on augmented system
    n = len(cols)
    rows = []  # each: (image dict, combo dict)
    piv = {}
    kern = []
    for i, img in enumerate(cols):
        img = {k: Fraction(v) for k, v in img.items() if v}
        comb = {i: Fraction(1)}
        while img:
            k = min(img)
            if k in piv:
                pimg, pcomb = piv[k]
                c = img[k] / pimg[k]
                for kk, vv in pimg.items():
                    nv = img.get(kk, 0) - c * vv
                    if nv: img[kk] = nv
                    else: img.pop(kk, None)
                for kk, vv in pcomb.items():
                    nv = comb.get(kk, 0) - c * vv
                    if nv: comb[kk] = nv
                    else: comb.pop(kk, None)
            else:
                piv[k] = (img, comb)
                break
        if not img:
            kern.append(comb)
    return kern

XI = ['z1', 'z2', 'w1', 'w2', 'w3']

def bulk_basis_A2(W):
    """W = (Wz1,Wz2,Ww1,Ww2,Ww3). elements (m, a, S): z^m w^(-a-1) xi^S, S sorted tuple of indices."""
    out = []
    for r in range(6):
        for S in combinations(range(5), r):
            s = [1 if i in S else 0 for i in range(5)]
            m = (W[0] + s[0], W[1] + s[1])
            a = tuple(-W[2 + I] - 1 - s[2 + I] for I in range(3))
            if min(m) < 0 or min(a) < 0: continue
            out.append((m, a, S))
    return out

def bulk_basis_A0(W):
    out = []
    for r in range(6):
        for S in combinations(range(5), r):
            s = [1 if i in S else 0 for i in range(5)]
            m = (W[0] + s[0], W[1] + s[1])
            b = tuple(W[2 + I] + s[2 + I] for I in range(3))
            if min(m) < 0 or min(b) < 0: continue
            out.append((m, b, S))
    return out

def contract(S, x):
    """d/dxi_x acting from the left on xi^S (ordered increasing). returns (sign, S') or None."""
    if x not in S: return None
    pos = S.index(x)
    return ((-1) ** pos, tuple(i for i in S if i != x))

def div_A2(el):
    m, a, S = el
    res = {}
    for x in S:
        sg, S2 = contract(S, x)
        if x < 2:
            if m[x] == 0: continue
            c = m[x]; m2 = list(m); m2[x] -= 1; key = (tuple(m2), a, S2)
        else:
            I = x - 2
            c = -(a[I] + 1); a2 = list(a); a2[I] += 1; key = (m, tuple(a2), S2)
        res[key] = res.get(key, 0) + sg * c
    return {k: v for k, v in res.items() if v}

def div_A0(el):
    m, b, S = el
    res = {}
    for x in S:
        sg, S2 = contract(S, x)
        if x < 2:
            if m[x] == 0: continue
            c = m[x]; m2 = list(m); m2[x] -= 1; key = (tuple(m2), b, S2)
        else:
            I = x - 2
            if b[I] == 0: continue
            c = b[I]; b2 = list(b); b2[I] -= 1; key = (m, tuple(b2), S2)
        res[key] = res.get(key, 0) + sg * c
    return {k: v for k, v in res.items() if v}

def F_times(el):
    """f xi_z1 xi_z2 * (z^m w^b xi^S) -> A2 element; zero unless b = 0 and S has no z."""
    m, b, S = el
    if any(b) or 0 in S or 1 in S: return {}
    # xi_z1 xi_z2 xi^S already ordered (0,1,...)
    return {(m, (0, 0, 0), (0, 1) + S): 1}

def bulk_coker(W):
    A2 = bulk_basis_A2(W)
    idx2 = {e: i for i, e in enumerate(A2)}
    imgs = [div_A2(e) for e in A2]
    ker2 = len(kernel_basis(imgs, None))
    WF = tuple(W[i] + 1 for i in range(5))
    A0 = bulk_basis_A0(WF)
    imgs0 = [div_A0(e) for e in A0]
    k0 = kernel_basis(imgs0, None)
    # F . ker -> A2 vectors -> del
    vecs = []
    for comb in k0:
        v = {}
        for i, c in comb.items():
            for k, cc in F_times(A0[i]).items():
                v[k] = v.get(k, 0) + c * cc
        v = {k: c for k, c in v.items() if c}
        if not v: continue
        dv = {}
        for k, c in v.items():
            for kk, cc in div_A2(k).items():
                dv[kk] = dv.get(kk, 0) + c * cc
        dv = {k: c for k, c in dv.items() if c}
        vecs.append(dv)
    return ker2 - rank(vecs), ker2, len(A2), len(A0)

# boundary de Rham complex of C^{2|3}
def bdy_basis(beta, p):
    """forms dz^Sz dtheta^a theta^J z^m with weight beta=(n1,n2,t1,t2,t3) and form degree p."""
    out = []
    for rz in range(3):
        for Sz in combinations(range(2), rz):
            m = tuple(beta[i] - (1 if i in Sz else 0) for i in range(2))
            if min(m) < 0: continue
            pa = p - rz
            if pa < 0: continue
            for J in product((0, 1), repeat=3):
                a = tuple(beta[2 + I] - J[I] for I in range(3))
                if min(a) < 0 or sum(a) != pa: continue
                out.append((Sz, a, J, m))
    return out

def d_bdy(el):
    Sz, a, J, m = el
    res = {}
    # dz_k d/dz_k : multiply dz_k on the left of dz^Sz (odd), derivative on z^m (even)
    for k in range(2):
        if m[k] == 0 or k in Sz: continue
        S2 = tuple(sorted(Sz + (k,)))
        sg = (-1) ** S2.index(k)
        m2 = list(m); m2[k] -= 1
        key = (S2, a, J, tuple(m2))
        res[key] = res.get(key, 0) + sg * m[k]
    # dtheta_I d/dtheta_I : d/dtheta_I passes dz^Sz (odd, |Sz|) and dtheta (even); acts on theta^J
    for I in range(3):
        if J[I] == 0: continue
        sg = (-1) ** len(Sz) * (-1) ** sum(J[:I])
        J2 = list(J); J2[I] = 0
        a2 = list(a); a2[I] += 1
        key = (Sz, tuple(a2), tuple(J2), m)
        res[key] = res.get(key, 0) + sg
    return {k: v for k, v in res.items() if v}

def bdy_B(beta, p):
    """dim of exact p-forms at weight beta = rank d: Omega^{p-1} -> Omega^p."""
    if p == 0: return 0
    return rank([d_bdy(e) for e in bdy_basis(beta, p - 1)])

def check_d2(beta):
    for p in range(0, sum(beta) + 3):
        for e in bdy_basis(beta, p):
            dd = {}
            for k, c in d_bdy(e).items():
                for kk, cc in d_bdy(k).items():
                    dd[kk] = dd.get(kk, 0) + c * cc
            if any(dd.values()): return False
    return True

if __name__ == "__main__":
    maxn = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    maxt = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    bad = 0; tot = 0; d2ok = True
    for n1 in range(maxn + 1):
        for n2 in range(n1, maxn + 1):
            for th in product(range(maxt + 1), repeat=3):
                if list(th) != sorted(th, reverse=True): continue
                beta = (n1, n2) + th
                d2ok &= check_d2(beta)
                W = (n1 - 1, n2 - 1) + tuple(-t - 1 for t in th)
                cok, ker2, dA2, dA0 = bulk_coker(W)
                maxp = sum(beta) + 3
                B = [bdy_B(beta, p) for p in range(maxp + 1)]
                pred = sum(B[2:]) + (1 if beta == (0, 0, 0, 0, 0) else 0)
                tot += 1
                flag = '' if cok == pred else '  MISMATCH'
                bad += cok != pred
                print(f"beta={beta}: bulk coker={cok} (ker del|A2={ker2}, dimA2={dA2}, dimA0={dA0}); "
                      f"B^p={B[1:]}  B^1={B[1]}  pred B^>=2+[0]={pred}{flag}")
    print(f"d^2=0 on all boundary bases: {d2ok};  slots {tot}, mismatches {bad}")
