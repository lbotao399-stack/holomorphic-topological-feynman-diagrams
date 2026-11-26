"""Momentum-space triangle rule (d = 2q - 2eps) on the prism: both cuts (rung = LA, triangle
edge = BE) reach a two-loop kite whose central line has power n + eps (FORCER master t1star05)
with an IR-divergent integer limit; the reduction does not close."""
import prop

PRISM = [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (2, 3), (2, 5), (3, 4), (4, 5)]
for cut in (0, 1):
    try:
        prop.cut_period(PRISM, [1, 1, 2, 1, 3, 4, 4, 4, 4], 6, cut=cut)
        print('cut', cut, 'reduced')
    except prop.NotReducible as e:
        print('cut', cut, 'NotReducible:', str(e)[:160])
print('star kites met:', prop.STAR_KITES)
