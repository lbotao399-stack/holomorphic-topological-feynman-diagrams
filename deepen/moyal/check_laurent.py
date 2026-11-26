"""Check the first Laurent coefficients in the full formal-N trace cochain algebra."""
from fractions import Fraction as F
from pathlib import Path
import sys
import pickle

DEEPEN = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(DEEPEN / "deform/code"))
import ce


def theta_reference():
    with (DEEPEN / "o3/code/data/theta.pkl").open("rb") as stream:
        return pickle.load(stream)


def cubic(degree, tensor):
    traces = {
        'f': {((0,1,2),):1,((0,2,1),):-1},
        'd': {((0,1,2),):1,((0,2,1),):1},
        's1': {((0,),(1,2)):1,((0,2),(1,)):1,((0,1),(2,)):1},
        's0': {((0,),(1,),(2,)):1},
    }
    bracket = ce.bracket_poly(0,1,3)
    polynomial = {(0,)*len(next(iter(bracket))): F(1)}
    for _ in range(degree):
        polynomial = ce.poly_mul(polynomial,bracket)
    colour = {word:{0:F(sign)} for word,sign in traces[tensor].items()}
    return ce.canon_raw(ce.raw_from(polynomial,colour),3)


def mark(left,right):
    return ce.bracket(left,3,right,3)


def laurent_terms():
    # All quantum vertices here have the common physical factor epsilon*pi*nu removed.
    d,s1,s0 = (cubic(1,t) for t in ('d','s1','s0'))
    projected = ce.lin((1,d),({-1:-2},s1),({-2:4},s0))
    minus_one = ce.lin(({2:F(-1,4)},projected))
    finite = ce.lin(({2:F(-3,20),0:F(-1,2)},cubic(2,'f')))
    linear = ce.lin(({2:F(-21,480),0:F(-1,12)},cubic(3,'d')),
                    ({1:F(3,80)},cubic(3,'s1')))
    return minus_one,finite,linear


def check():
    minus_one,finite,linear = laurent_terms()
    classical = [cubic(0,'f'),ce.scale(cubic(1,'d'),F(1,2)),ce.scale(cubic(2,'f'),F(1,8))]
    data = theta_reference()
    bare = [data['op'],ce.scale(data['th']['d'],F(1,2))]
    assert not ce.d_ce(minus_one,3)
    assert not ce.lin((-1,ce.d_ce(finite,3)),(1,mark(classical[1],minus_one)),(1,bare[0]))
    assert not ce.lin((-1,ce.d_ce(linear,3)),(1,mark(classical[1],finite)),
                      (1,mark(classical[2],minus_one)),(1,bare[1]))
    assert not ce.lin((1,mark(classical[1],classical[1])),(-2,ce.d_ce(classical[2],3)))
    assert not ce.lin((1,mark(classical[0],finite)),(1,ce.d_ce(finite,3)))
    print("e^-1, e^0, e^1 equations passed at formal N; Moyal pole redefinitions checked.")


if __name__=='__main__':
    check()
