"""H^3 of g[[z1,z2]] mod translations (invariant cochains, trace basis, formal N) at the weights
(1,1), (2,0), (0,2): C^2 = 0, Z^3 = ker(d: C^3 -> C^4), exact; GL(2) part [p1p2] T.
Saves Z^3 bases to h3_data.pkl."""
import sys, time, pickle
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/deform/code')
import ce, la, spaces


def show(c, k):
    terms = []
    for (t, m), v in sorted(c.items(), key=lambda x: (x[0][0], x[0][1])):
        mon = "".join(("p%d%s" % (s + 1, "ab"[cc])) + ("^%d" % m[2 * s + cc] if m[2 * s + cc] > 1 else "")
                      for s in range(k - 1) for cc in range(2) if m[2 * s + cc])
        terms.append("%s %s %s" % (v if not isinstance(v, dict) else v, mon or "1", ce.tau_str(t)))
    return " + ".join(terms)


t0 = time.time()
DATA = {}
for sl in (False, True):
    name = "sl_N" if sl else "gl_N"
    for deg in [(1, 1), (2, 0), (0, 2)]:
        C2, _ = spaces.basis(2, deg, sl)
        C3, keys3 = spaces.basis(3, deg, sl)
        keys4 = spaces.keys_of(4, deg, sl)
        dC3 = [spaces.const_part(ce.d_ce(spaces.to_np(c), 3, sl)) for c in C3]
        _, M = la.matrix_from(dC3, keys4)
        cert = la.certified_rank(M)
        Z3 = []
        for v in cert['kernel']:
            z = {}
            for j, cj in enumerate(v):
                if cj:
                    for k, x in C3[j].items():
                        z[k] = z.get(k, 0) + cj * x
            Z3.append({k: x for k, x in z.items() if x})
        closed = all(ce.d_ce(spaces.to_np(z), 3, sl) == {} for z in Z3)
        print("%s deg %s: dim C^2 = %d, dim C^3 = %d, rank d3 = %d (det %s), dim Z^3 = dim H^3 = %d (closed exactly: %s)" % (
            name, deg, len(C2), len(C3), cert['rank'], cert['det'], len(Z3), closed), flush=True)
        DATA[(name, deg)] = dict(C3=C3, Z3=Z3, keys3=keys3)
    # GL(2)-equivariant part at (1,1): [p1p2] T with T in Sym^3 invariants
    u = ce.bracket_poly(0, 1, 3)
    Ts = {"d=tr(1{23})": {((0, 1, 2),): {0: Fr(1)}, ((0, 2, 1),): {0: Fr(1)}},
          "s1=tr1tr23+tr2tr13+tr3tr12": {((0,), (1, 2)): {0: Fr(1)}, ((0, 2), (1,)): {0: Fr(1)}, ((0, 1), (2,)): {0: Fr(1)}},
          "s0=tr1tr2tr3": {((0,), (1,), (2,)): {0: Fr(1)}}}
    if sl:
        Ts = {"d=tr(1{23})": Ts["d=tr(1{23})"]}
    Js = {}
    for nm, T in Ts.items():
        J = ce.canon_raw(ce.raw_from(u, T), 3, sl=sl)
        J2 = ce.scale(ce.alt_raw(ce.raw_from(u, T), 3, sl=sl), Fr(1, 6))
        same = ce.clean(ce.lin((1, J), (-1, J2))) == {}
        print("  %s: J = [p1p2] T alternating: %s, alt([p1p2]T) = 6 J: %s, d J = 0: %s" % (
            nm, ce.is_alternating(J, 3, sl), same, ce.d_ce(J, 3, sl) == {}))
        Js[nm] = spaces.const_part(J)
    # non-symmetric T give alt([p1p2] T) = 6 [p1p2] Sym(T): check on all 6 taus
    for tau in spaces.taus(3, sl):
        T = {tau: {0: Fr(1)}}
        A = ce.alt_raw(ce.raw_from(u, T), 3, sl=sl)
        sym = {}
        for s, sg in ce.perms(3):
            ce.add_term(sym, ce.relabel(tau, {i: (s[i],) for i in range(3)}), {0: Fr(1)})
        B = ce.canon_raw(ce.raw_from(u, sym), 3, sl=sl)
        print("    alt([p1p2] %s) = [p1p2] * sum_s s.T : %s" % (ce.tau_str(tau), ce.clean(ce.lin((1, A), (-1, B))) == {}))
    # GL(2) part inside Z^3(1,1)
    Z = DATA[(name, (1, 1))]['Z3']
    keys3 = DATA[(name, (1, 1))]['keys3']
    _, MZ = la.matrix_from(Z, keys3)
    cz = la.certified_rank(MZ)
    for nm, J in Js.items():
        print("  %s in Z^3(1,1): %s" % (nm, la.in_span(MZ, la.column(J, keys3), cz)[0]))
    DATA[(name, 'GL2')] = Js
    print("  Z^3(1,1) basis:")
    for i, z in enumerate(Z):
        print("   z%d = %s" % (i, show(z, 3)))
pickle.dump(DATA, open((__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/deform/code/h3_data.pkl', 'wb'))
print("total %.0fs" % (time.time() - t0))
