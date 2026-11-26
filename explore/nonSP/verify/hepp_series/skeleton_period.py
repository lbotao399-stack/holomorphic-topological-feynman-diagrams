"""Exact Gegenbauer expansion of a reconstructed four-loop nine-edge numerator."""
import argparse
import json
from fractions import Fraction as F
from pathlib import Path

from period import HeppSeries


def load(path):
    data = json.loads(path.read_text())
    poly = {tuple(a+b for a,b in zip(e, data['h'])): F(c, data['denominator'])
            for e,c in data['terms']}
    return HeppSeries(data['edges'], F(11,2), poly)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('numerator',type=Path)
    parser.add_argument('order',type=int)
    opts = parser.parse_args()
    engine = load(opts.numerator)
    print('cotrees',engine.tree_count,'monomials',len(engine.numerator),flush=True)
    signs={int(c>0)-int(c<0) for c in engine.numerator.values()}
    if len(signs)==1 and opts.order>=3:
        from moment_bounds import enclose
        sign=signs.pop()
        moments=[sign*m for k,m,b in engine.moments(3,permutations='matroid')]
        lo,hi=enclose(moments,engine.alpha,engine.tree_count)
        lower,upper=(lo,hi) if sign>0 else (-hi,-lo)
        print('Hermite enclosure with integer endpoints:',lower//1,-((-upper)//1),flush=True)
    for k,t,p,b in engine.terms(opts.order,permutations='graph'):
        print(k,'term',t,'partial',p,'tail <=',b,flush=True)
