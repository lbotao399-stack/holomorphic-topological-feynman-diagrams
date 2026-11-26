"""All-loop local matrix element of the ordered six-dimensional quantum differential."""
from collections import Counter
from fractions import Fraction
from itertools import permutations
from math import factorial, prod
from pathlib import Path
from random import Random
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'code'))
from htk import Graph, brute
from rooted_trees import graphs
from wick_weights import wick_sign


def path(n):
    root = 2*n+1
    edges = []
    for j in range(n):
        edges.extend([(root,2*j),(2*j,2*j+1),(2*j+1,2*j+2)])
    return edges+[(2*n,root)]


def period(n):
    return Fraction(2**(n+1),factorial(n+2)*factorial(n)**2*factorial(n+1)**2)


def original_times(ratios):
    effective = jacobian = 1
    tail = []
    for w,s1,s2 in reversed(ratios):
        jacobian *= effective**3/(w**4*(1-w)**2)
        series, direct = effective/w, effective/(1-w)
        tail = [series*(s1-s2),series*s2,direct]+tail
        effective = series*(1-s1)
    return [effective]+tail, jacobian


def check_density(n,rng):
    graph = Graph(2*n+2,path(n))
    sources = {0:(1,0,0),2*n:(0,0,1)}
    sources.update({2*j+1:(0,1,0) for j in range(n)})
    ratios = []
    for _ in range(n):
        s1,s2 = sorted(rng.sample(range(1,16),2), reverse=True)
        ratios.append(tuple(Fraction(k,16) for k in (rng.randrange(1,16),s1,s2)))
    times, jacobian = original_times(ratios)
    raw = graph.sigma()*brute.jacobian_density(graph,times,sources)
    expected = prod((1-w)**2*w**j*(1-s1)**(n-j-1) for j,(w,s1,_) in enumerate(ratios))
    assert raw*jacobian == expected
    integrated = prod(Fraction(2*factorial(j),factorial(j+3)*(n-j)*(n-j+1)) for j in range(n))
    assert integrated == period(n)
    assert wick_sign(path(n),graph.root)*(-1)**(n-1) == (-1)**(n*(n-1)//2)
    # Impose total output momentum zero by changing the last middle momentum.
    sources[2*n-1] = (-1,1-n,-1)
    w,s1,s2 = ratios[-1]
    first, middle = 1,0
    for _,t1,t2 in ratios[:-1]:
        first, middle = (1-t1)*first,(1-t1)*middle+1-t2
    correction = (-s2*(1-w)/w)**(n-1)*(Fraction(middle,first)-(n-1))
    assert graph.sigma()*brute.jacobian_density(graph,times,sources)*jacobian == expected*correction
    mean = 0
    previous = 1
    for k in reversed(range(1,n)):
        mean += previous*Fraction(k+1,k)
        previous *= Fraction(k+2,k)
    assert mean-(n-1) == Fraction((n-1)*(3*n-2),4)


def parity(order):
    return (-1)**sum(x>y for j,x in enumerate(order) for y in order[j+1:])


def bracket(x,y):
    """Sparse matrix-unit commutator, with no dimension or color identities imposed."""
    z = Counter()
    for (i,j),v in x.items():
        for (k,l),w in y.items():
            if j == k:
                z[i,l] += v*w
            if l == i:
                z[k,j] -= w*v
    return {k:v for k,v in z.items() if v}


def color_function(edges,root):
    spokes = [e for e,edge in enumerate(edges) if root in edge]
    outputs = [v for v in range(root) if sum(v in edge for edge in edges)==2]
    legs = []
    for v in range(root):
        inc, out = [], []
        for e,(tail,head) in enumerate(edges):
            if v not in (tail,head):
                continue
            other = head if tail == v else tail
            child = ('leaf',spokes.index(e)) if other == root else other
            (inc if head==v else out).append(child)
        legs.append(inc+([('leaf',len(spokes)+outputs.index(v))] if v in outputs else [])+out)
    anchor = next(v for v in edges[spokes[0]] if v != root)
    def rooted(v, parent):
        if isinstance(v,tuple):
            return v[1]
        i = legs[v].index(parent)
        return (rooted(legs[v][(i+1)%3],v),rooted(legs[v][(i+2)%3],v))
    tree = rooted(anchor,('leaf',0))
    def evaluate(leaves):
        def descend(branch):
            if isinstance(branch,int):
                return leaves[branch]
            return bracket(descend(branch[0]),descend(branch[1]))
        result = descend(tree)
        return sum(v*result.get((j,i),0) for (i,j),v in leaves[0].items())
    return evaluate


def check_colors(n):
    m = n+1
    units = [{(j,(j+1)%(2*m+1)):1} for j in range(2*m+1)]
    contributing = []
    for edges in graphs(n):
        evaluate = color_function(edges,2*n+1)
        entries = 0
        for ins in permutations(range(m)):
            for outs in permutations(range(m,2*m+1)):
                entries += bool(evaluate([units[j] for j in ins+outs]))
        if entries:
            contributing.append(entries)
    assert contributing == [2]
    ordered = color_function(path(n),2*n+1)
    forward = tuple(reversed(range(m)))+tuple(range(m,2*m+1))
    reverse = tuple(range(m))+tuple(reversed(range(m,2*m+1)))
    coefficient = 0
    for order,reflect in ((forward,False),(reverse,True)):
        color = ordered([units[j] for j in order])
        coefficient += color*parity(order[:m])*parity(order[m:])*((-1)**n if reflect else 1)
    derivative = (-1)**(n*(n+1)//2)
    weight = Fraction(wick_sign(path(n),2*n+1),2)
    assert coefficient*derivative*weight == (-1)**(n*(n-1)//2)


if __name__ == '__main__':
    rng = Random(71)
    for n in range(1,11):
        check_density(n,rng)
        if n<=3:
            check_colors(n)
        local = (-1)**(n*(n-1)//2)*period(n)
        integrated = local*Fraction((-1)**(n-1)*(n-1)*(3*n-2),4)
        print(n, 'local and zero-total-momentum coefficients in units i/(8 pi^3)^n', local, integrated)
