"""Compare the coefficient detector against the project's CE differential with matrix-unit inputs."""
from math import prod
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'deform' / 'code'))
import ce

p = [(1,0), (0,1), (1,1), (-2,-2)]
coords = sum((list(a) for a in p[:-1]), [])
for w in range(2, 11):
    poly = {tuple([0] * 6): 1}
    for _ in range(w):
        poly = ce.poly_mul(poly, ce.bracket_poly(0, 1, 3))
    color = {((0,1,2),): {0: 1}, ((0,2,1),): {0: (-1)**(w+1)}}
    cubic = ce.canon_raw(ce.raw_from(poly, color), 3, sl=True)
    quartic = ce.d_ce(cubic, 3, sl=True)
    value = sum(coeff[0]*prod(a**i for a, i in zip(coords, m))
                for (word, m), coeff in quartic.items() if word == ((0,1,2,3),))
    assert value == 2**w + (-1)**w - 1
    print(w, value)
