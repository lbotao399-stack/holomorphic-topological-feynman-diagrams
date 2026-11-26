"""Completion identity: Pi(Gc) = Gamma(q) I_q(Gc - c; q - b) / prod Gamma(b_e) independent of c.
For every convergent orbit of a tabulated core G (K4, W4, W5, prism, K33) at level q: complete,
evaluate every decompletion Gc - c with all b_e in [1, q-1] that cores.period can evaluate, compare.
Also: cut criterion (cut(W) > 2q, 2 <= |W| <= V-2) vs gutil.convergent on every decompletion."""
import sys
import time
from collections import Counter

import engine as en
import gutil

K4 = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
W4 = [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (2, 3), (3, 4), (4, 1)]
W5 = [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (1, 2), (2, 3), (3, 4), (4, 5), (5, 1)]
PRISM = list(en.cores.PRISM)
K33 = list(en.cores.K33)
CASES = [('K4', K4, (3, 4, 5, 6, 7)), ('W4', W4, (3, 4, 5, 6)), ('W5', W5, (3, 4)),
         ('prism', PRISM, (3, 4, 5)), ('K33', K33, (3, 4, 5))]


def main(CASES):
  for name, G, qs in CASES:
      n = gutil.cycle_rank(G)
      autos = gutil.edge_automorphisms(tuple(G))
      for q in qs:
          t = time.time()
          reps = set()
          for c in gutil.compositions(n * q, len(G), 1, q - 1):
              if gutil.convergent(tuple(G), list(c), q):
                  reps.add(gutil.orbit_rep(c, autos))
          pairs = bad = graphs = 0
          conv_bad = conv_n = 0
          types = Counter()
          for al in sorted(reps):
              V, E, inf = en.complete(G, al, q)
              assert en.balanced(V, E, q)
              cutc = en.cut_convergent(V, E, q)
              vals = []
              for c in range(V):
                  H = en.decomplete(V, E, c, q)
                  if H is None:
                      continue
                  conv_n += 1
                  if gutil.convergent(tuple(H[0]), H[1], q) != cutc:
                      conv_bad += 1
                  r = en.pi_from_H(H, q)
                  if r is not None:
                      vals.append((c, r[0], r[1]))
              if len(vals) >= 2:
                  graphs += 1
                  for c, v, ty in vals[1:]:
                      pairs += 1
                      types[(vals[0][2], ty)] += 1
                      if v != vals[0][1]:
                          bad += 1
                          print('MISMATCH', name, q, al, vals[0], (c, v, ty), flush=True)
          print('%-5s q=%d orbits %d  graphs with >=2 tabulated decompletions %d  pairs %d  mismatches %d'
                '  | cut criterion vs gutil.convergent: %d decompletions, %d disagreements  (%.1fs)'
                % (name, q, len(reps), graphs, pairs, bad, conv_n, conv_bad, time.time() - t), flush=True)
          print('      pair types:', dict(types), flush=True)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        CASES = [c for c in CASES if c[0] in sys.argv[1:]]
    main(CASES)
