"""Exact signed maximal minors for all eight E6 benchmark topologies."""
from fractions import Fraction as F
from pathlib import Path
from random import Random
import sys
import sympy as s

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from graphs import Graph, K4SUB, K33, PRISM
from htk import brute


def chamber(graph, times, source):
    current = s.Matrix(brute.currents(graph,times,source))
    C = s.Matrix(graph.cycle_matrix())
    psi = s.Matrix.vstack(*(C*s.diag(*current.col(i)) for i in range(current.cols)))
    E = graph.ne
    minors = s.Matrix([(-1)**e*psi[:,[i for i in range(E) if i!=e]].det() for e in range(E)])
    assert psi*minors == s.zeros(E-1,1)
    if not minors[0]:
        assert minors == s.zeros(E,1)
        return 0
    assert minors/minors[0] == s.Matrix(times)/times[0]
    U = (C*s.diag(*times)*C.T).det()
    density = F(minors[0]/(times[0]*U**current.cols))
    assert density == brute.density(graph,times,source).rational()
    return density


if __name__=='__main__':
    rng = Random(32864)
    for name,edges,h in [(n,e,3) for n,e in K4SUB]+[('prism',PRISM,2),('K33',K33,2)]:
        graph = Graph(max(max(e) for e in edges)+1,edges)
        for _ in range(2):
            examples = 0
            while not examples:
                source = [[rng.randint(-3,3) for _ in range(h)] for _ in range(graph.root)]
                times = [F(rng.randint(1,9)) for _ in edges]
                examples = chamber(graph,times,source)
        print(name,'two exact projective inverses and density signs')
    graph = Graph(5,[(0,1),(0,2),(0,3),(1,4),(4,2),(1,3),(2,3)])
    source = [(1,0),(0,0),(0,1),(-1,0)]
    minus = [1,1,1,1,1,1,2]
    plus = [1,1,1,1,1,2,1]
    assert chamber(graph,minus,source) == -F(6,52521875)
    assert chamber(graph,plus,source) == F(6,52521875)
    print('Paired sources occupy both signed chambers.')
