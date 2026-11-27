"""ESW holomorphic (T = 0) twists as quadratic Lie superalgebras g_tot on C^H  [src ESW arXiv:2002.10517 Part II;
local TeX conversion 09_logs/taxonomy_notion_import/generated/*.md].

  L_R  = g |x PiR                         (Map(C^H, R/g):  c in g, gamma in R even letters)
  L_U  = g + PiU + g^v                    ((U (x) K^{1/2})//g:  c, phi in U, Koszul generator xi in g^v)
  T[1] = (x) Lambda[eps]                  (odd tangent shift);  C_Dol, g/g: further odd generators
  T*_p = BF double with pairing parity p = H + 1 mod 2

  d   N        ESW thm      H  g_tot
  10  (1,0)    4.4          5  g                                    (hCS, CS class)
  8   1        6.1          4  T*_1 g                = g (x) Lambda[th]
  6   (1,1)    8.3          3  T*_0 (g (x) Lambda[e]) = g (x) Lambda[th1,th2]  (Berezin form = 6d note)
  6   (1,0)    8.1, 8.2     3  L_U                                  (U = 0: T*_0 g)
  4   1        10.1         2  T*_1 L_R
  4   2        10.2         2  T*_1 L_U
  4   4        10.9         2  T*_1 (g (x) Lambda[e1,e2]) = g (x) Lambda[th1,th2,th3] (Berezin form = LLC)
  2   (2,2)    12.1         1  T*_0 (L_R (x) Lambda[e])
  2   (0,2)    12.16        1  T*_0 L_R
  2   (4,4)    12.7         1  T*_0 (L_U (x) Lambda[e])
  2   (8,8)    12.3.1       1  T*_0 (g (x) Lambda[e1,e2,e3]) = g (x) Lambda[th1..th4]
"""
import superalg as sa


def gauge(name):
    if name == 'su2':
        return sa.su2_eps()
    if name.startswith('sl'):
        return sa.sl(int(name[2:]))
    if name.startswith('so'):
        return sa.so(int(name[2:]))
    if name.startswith('sp'):
        return sa.sp(int(name[2:]))
    if name.startswith('gl'):
        return sa.gl(int(name[2:]))
    raise ValueError(name)


def rep(g, spec):
    """spec: 'fund', 'adj', 'fund*k' (k copies), 'fund+adj', ..."""
    parts = []
    for item in spec.split('+'):
        k = 1
        if '*' in item:
            item, k = item.split('*')
            k = int(k)
        r = sa.fundamental(g) if item == 'fund' else sa.adjoint_rep(g) if item == 'adj' else \
            sa.dual_rep(sa.fundamental(g)) if item == 'afund' else None
        if r is None:
            raise ValueError(item)
        parts += [r] * k
    return sa.direct_sum(*parts) if len(parts) > 1 else parts[0]


def hyper(g, spec):
    """spec 'T*R' style: U = R + R^*; 'half:fund' : U = fundamental of sp(2N) with Om = J."""
    if spec.startswith('half:'):
        assert hasattr(g, 'J') and spec == 'half:fund'
        return sa.fundamental(g), g.J.copy()
    return sa.cotangent_hyper(rep(g, spec))


def build(theory, gname, matter=None):
    g = gauge(gname)
    E = sa.lie_even
    if theory == '10d':
        return 5, E(g)
    if theory == '8d':
        return 4, sa.double(E(g, False), 1)
    if theory == '6d11':
        return 3, sa.berezin_tensor(g, 2)
    if theory == '6d11b':
        return 3, sa.double(sa.tensor_lambda(E(g, False), 1), 0)
    if theory == '6d10':
        if not matter:
            return 3, sa.double(E(g, False), 0)
        rho, Om = hyper(g, matter)
        return 3, sa.hyper_LU(g, rho, Om)
    if theory == '4dN1':
        if not matter:
            return 2, sa.double(E(g, False), 1)
        return 2, sa.double(sa.semidirect(g, rep(g, matter), 'r'), 1)
    if theory == '4dN2':
        if not matter:
            L = sa.double(E(g, False), 0)          # L_U with U = 0
            L.w = None
            return 2, sa.double(L, 1, suffix='*')
        rho, Om = hyper(g, matter)
        return 2, sa.double(sa.hyper_LU(g, rho, Om), 1, suffix='*')
    if theory == '4dN4':
        return 2, sa.berezin_tensor(g, 3)
    if theory == '4dN4b':
        return 2, sa.double(sa.tensor_lambda(E(g, False), 2), 1)
    if theory == '2d22':
        return 1, sa.double(sa.tensor_lambda(sa.semidirect(g, rep(g, matter), 'r'), 1), 0)
    if theory == '2d02':
        return 1, sa.double(sa.semidirect(g, rep(g, matter), 'r'), 0)
    if theory == '2d44':
        rho, Om = hyper(g, matter)
        return 1, sa.double(sa.tensor_lambda(sa.hyper_LU(g, rho, Om), 1), 0)
    if theory == '2d88':
        return 1, sa.berezin_tensor(g, 4)
    raise ValueError(theory)


def catalogue_small():
    specs = [('10d', 'su2', None), ('8d', 'su2', None), ('6d11', 'su2', None), ('6d11b', 'su2', None),
             ('6d10', 'sl2', None), ('6d10', 'sl2', 'fund'), ('4dN1', 'sl2', None), ('4dN1', 'sl2', 'fund'),
             ('4dN2', 'sl2', None), ('4dN2', 'sl2', 'fund'), ('4dN4', 'su2', None), ('4dN4b', 'su2', None),
             ('2d22', 'sl2', 'fund'), ('2d02', 'sl2', 'fund'), ('2d44', 'sl2', 'fund'), ('2d88', 'su2', None),
             ('6d10', 'sp2', 'half:fund')]
    out = {}
    for th, g, m in specs:
        H, alg = build(th, g, m)
        out['%s/%s/%s' % (th, g, m)] = dict(H=H, alg=alg)
    return out
