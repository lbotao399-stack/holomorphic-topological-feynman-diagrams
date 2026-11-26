"""Extra validation of the character formula (t7_sizes.dimC) against explicit enumeration for n = 4, w = 1,
G = -5..-3, all weights.  Writes check_char_n4.log."""
import time
from cohom import Space, all_weights
from t7_sizes import dimC
out = open('check_char_n4.log', 'w')
n = 4
S = Space(n, 1, 2)
ok = True; cnt = 0; tot = 0
for G in range(-5, -2):
    for t in all_weights(S, G):
        sl = S.slot(t)
        for N in range(0, 30):
            d = len(sl.get(N, []))
            dd = dimC(n, N, 'g', t) + (dimC(n, N - 1, 'b', t) if N >= 1 else 0)
            if d or dd:
                cnt += 1; tot += d
                if d != dd:
                    ok = False
                    out.write('MISMATCH %s %s %d %d %d\n' % (G, t, N, d, dd))
        S.drop(t)
s = 'n=4 w=1 G=-5..-3: %d nonzero slots (total dim %d): character formula == enumeration: %s' % (cnt, tot, ok)
print(s); out.write(s + '\n')
