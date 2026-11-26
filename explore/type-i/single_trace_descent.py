"""Solve the single-trace projected Ward equation for SO(N), including the mixed hexagon tensor."""
from fractions import Fraction as F
import json
import sys
import sympy as sp

from open_ce import coordinates,differential
from open_theta import HERE
from mixed_hexagon import single_trace


def main():
    pure=[{} for _ in range(3)]
    for row in json.loads((HERE/'open_theta.json').read_text())['sum']:
        if row['lengths']==[8]:
            pure[row['rank']]={(tuple(a),tuple(b)):F(c) for a,b,c in row['terms']}
    inputs=[s for s in coordinates() if s[0]==(7,)]
    exact=[{q:v for (_,q),v in differential(*s).items()} for s in inputs]
    mixed=single_trace()
    rows=sorted(set().union(*[set(d) for d in [*pure,*exact,mixed]]))
    available=sp.Matrix([[d.get(q,0) for d in exact] for q in rows])
    _,ids=available.rref()
    images=available[:,list(ids)].row_join(sp.Matrix([mixed.get(q,0) for q in rows]))
    targets=sp.Matrix([[d.get(q,0) for d in pure] for q in rows])
    _,pivots=images.T.rref()
    assert len(ids)==3 and len(pivots)==4
    solution=images[list(pivots),:].inv()*targets[list(pivots),:]
    n=sp.Symbol('N')
    rank=sp.Matrix([1,n,n**2])
    residual=(targets-images*solution)*rank
    assert sp.Poly(sp.gcd_list([e for e in residual if e]),n).monic().as_expr()==n-32
    x=solution*sp.Matrix([1,32,1024])
    assert x[-1]==sp.Rational(1,8640)
    assert images*x==targets*sp.Matrix([1,32,1024])
    result={'orthogonal_rank':32,'mixed_tensor_subtraction':str(x[-1]),
            'pure_gauge_preimages':[[[list(s) for s in inputs[j]],str(x[i])] for i,j in enumerate(ids) if x[i]]}
    path=HERE/'single_trace_preimages.json'
    if '--save' in sys.argv:path.write_text(json.dumps(result,separators=(',',':'))+'\n')
    assert result==json.loads(path.read_text())
    print(result)


if __name__=='__main__':
    main()
