"""(c) I_2 for sl_N (formal N, kappa = tr), j = a D, order-hbar^2 equation  eps(-d i_2 + (1/2){j,j}) + o_2 = 0,
o_2 = pi nu o', c_O = pi nu/48, a^2 = -12 N^2 eps c_O = -eps pi nu N^2/4.
Write i_2 = eps pi nu iota; then  d iota = -(N^2/4)(1/2){D,D} + o'.
Checks (exact, formal N):
  (1/2){D,D} = F_223/(12N^2) + d eta_N,  eta_N = ((N^2+4)/(2N^2)) [p1p2]^2 f        (deform d7)
  F_133 + F_223 = 2 N^2 beta,  beta = d([p1p2]^2 f)                                   (beta_Killing/2 with Killing = 2N tr)
  iota = -((3N^2+10)/20) [p1p2]^2 f solves the equation;  d injective on C^3_(2,2) (uniqueness).
Also gl_N: the same identities checked where they hold.  Saves data/i2_sl.pkl."""
import time
from common import *

t0 = time.time()
th = load('theta.pkl')
obs = spaces.obstruction_cochains()
res = {}
for sl in (True, False):
    nm = 'sl_N' if sl else 'gl_N'
    F133 = ce.restrict_sl(obs[(1, 3, 3)]['F']) if sl else obs[(1, 3, 3)]['F']
    F223 = ce.restrict_sl(obs[(2, 2, 3)]['F']) if sl else obs[(2, 2, 3)]['F']
    op = ce.restrict_sl(th['op']) if sl else th['op']
    D = J_cochain('d', sl)
    u2f = f_bracket_sq(sl)
    beta = ce.d_ce(u2f, 3, sl)
    half = ce.scale(ce.bracket(D, 3, D, 3, sl=sl), Fr(1, 2))
    # eta_N = ((N^2+4)/(2N^2)) u2f  -> Npoly coefficient {0: 1/2, -2: 2}
    etaN = ce.lin(({0: Fr(1, 2), -2: Fr(2)}, u2f))
    lhs = ce.lin((1, half), ({-2: Fr(-1, 12)}, F223), (-1, ce.d_ce(etaN, 3, sl)))
    print("%s: (1/2){D,D} - F_223/(12N^2) - d eta_N == 0: %s" % (nm, is_zero(lhs)))
    s = ce.lin((1, F133), (1, F223), ({2: Fr(-2)}, beta))
    print("%s: F_133 + F_223 - 2N^2 d([p1p2]^2 f) == 0: %s" % (nm, is_zero(s)))
    if sl:
        iota = ce.lin(({2: Fr(-3, 20), 0: Fr(-1, 2)}, u2f))
        eq = ce.lin((1, ce.d_ce(iota, 3, sl)), ({2: Fr(1, 4)}, half), (-1, op))
        print("%s: d iota == -(N^2/4)(1/2){D,D} + o',  iota = -((3N^2+10)/20)[p1p2]^2 f: %s" % (nm, is_zero(eq)))
        C3, k3 = spaces.basis(3, (2, 2), sl)
        k4 = spaces.keys_of(4, (2, 2), sl)
        img = [spaces.const_part(ce.d_ce(spaces.to_np(c), 3, sl)) for c in C3]
        _, M = la.matrix_from(img, k4)
        cert = la.certified_rank(M)
        print("%s: dim C^3_(2,2) = %d, rank d = %d (minor det %s): d injective, iota unique: %s; dim C^2_(2,2) = %d"
              % (nm, len(C3), cert['rank'], cert['det'], cert['rank'] == len(C3), len(spaces.basis(2, (2, 2), sl)[0])))
        print("  [p1p2]^2 f =", show(u2f, 3))
        res['iota'] = iota
        res['u2f'] = u2f
save(res, 'i2_sl.pkl')
print("time %.1fs" % (time.time() - t0))
