"""Tree and marked two-loop parts of the sl_N order-four coefficient probe."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations, permutations, product
from math import prod
from pathlib import Path
import sys
import sympy as sp

root = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(root/'moyal'))
from theta_moyal import Theta, ribbon, ce, calculate
from check_laurent import cubic

P = [(1,0),(0,1),(1,1),(-2,-2)]
N = sp.Symbol('N')
word = ((0,1,2,3),)
terms = {'d':(0,1), 's1':(-1,-2), 's0':(-2,4)}


def value(cochain):
    coords = sum((list(p) for p in P[:-1]),[])
    return sp.expand(sum(prod(x**a for x,a in zip(coords,m))*sum(v*N**i for i,v in powers.items())
                         for (t,m),powers in cochain.items() if t == word))


def oriented_sum(colour, monomials):
    result = defaultdict(F)
    for s in permutations(range(4)):
        q = [P[i] for i in s]
        b = [q[i][0]*q[j][1]-q[j][0]*q[i][1] for i,j in ((0,1),(0,2),(1,2))]
        kernel = sum(coeff*prod(x**a for x,a in zip(b,m)) for m,coeff in monomials.items())
        for t,by_rank in colour.items():
            if ribbon.canon_tau([tuple(s[i] for i in w) for w in t]) == word:
                for rank,coeff in by_rank.items():
                    result[rank] += ribbon.perm_sign(s)*kernel*coeff
    return result


def marked_thetas(shape):
    theta = Theta(shape)
    original = theta.brackets[0]*(theta.brackets[1]+theta.brackets[2])
    for marks in [(i,i) for i in range(6)] + list(combinations(range(6),2)):
        i,j = marks
        poly = original*theta.vertices[i]*theta.vertices[j]
        integrated = defaultdict(F)
        for m,coeff in poly.items():
            integrated[m[:3]] += int(coeff)*theta.moment(m[3:])*theta.prefactor/4
        if i == j:
            colours = [({}, 0, 1, (0,))]
        else:
            colours = [({i:s,j:t},terms[s][0]+terms[t][0],terms[s][1]*terms[t][1],
                        (1,2) if s==t=='d' else (1,)) for s,t in product(terms,repeat=2)]
        for assignments,power,coeff,slots in colours:
            c = ribbon.colour_weight(6,theta.edges,theta.ext,assignments)
            yield c,integrated,power,coeff,slots


def theta_values(shape):
    output = [defaultdict(F) for _ in range(3)]
    for c,integrated,power,coeff,slots in marked_thetas(shape):
        for k,v in oriented_sum(c,integrated).items():
            for slot in slots:
                output[slot][k+power] += coeff*v
    return [sp.expand(sum(value*N**k for k,value in comp.items())) for comp in output]


def project():
    tree_dd = value(ce.bracket(cubic(1,'d'),3,cubic(3,'d'),3,sl=True))
    tree_ff = value(ce.bracket(cubic(2,'f'),3,cubic(2,'f'),3,sl=True))/2
    parts = [theta_values(s) for s in ((1,3,3),(2,2,3))]
    single, double, unprojected = [sum(row[j] for row in parts) for j in range(3)]
    moyal = calculate(2)
    x = (1,1,-1)  # [01], [02], [12] on P
    reference = sum(v*N**i*prod(a**b for a,b in zip(x,m))
                    for (t,i,m),v in moyal.items() if t == word)
    assert sp.expand(reference-single/8-unprojected/4)==0
    # Divide by (epsilon pi nu)^2. Then a^2/(epsilon pi nu)=-N^2/4.
    a2 = -N**2/4
    k2 = -(3*N**2+10)/20
    k3 = -(2*N**2+5)/30  # i3 coefficient divided by a epsilon pi nu
    tree = sp.expand(a2*k3*tree_dd+k2**2*tree_ff)
    theta = sp.expand(k2*single+a2*double)
    print('tree bases:',tree_dd, tree_ff)
    for shape, row in zip(('133','223'),parts):
        print(shape,'one f2, two projected d, two unprojected d:',*row)
    coefficient = sp.factor((tree+theta)/16)
    print('tree:',tree,'two loops:',theta,'coefficient / 16:',coefficient)
    return coefficient


if __name__=='__main__':
    project()
