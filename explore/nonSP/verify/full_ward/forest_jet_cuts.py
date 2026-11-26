"""Independent original-graph checks for the two zero-bare-kernel topologies."""
from fractions import Fraction
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'prism_forcer'))
from four_loop import evaluate
from catalogue import graph_catalogue
from forest_jets import evaluate as forest


def check(folder,executable):
    graphs=graph_catalogue()
    source=[(0,0),(0,1),(1,0),(-3,-2),(0,0),(2,1)]
    edge=graphs[6].index((0,2))
    ref=forest(graphs[6],2,3,source,[(edge,1)])
    assert ref==-Fraction(1,259200)
    assert evaluate(folder,executable,graphs[6],source,[(edge,1)])==ref
    source=[(0,0),(0,1),(1,0),(0,0),(-1,-1),(0,0)]
    edge=graphs[10].index((2,3))
    ref=forest(graphs[10],2,4,source,[(edge,0)])
    assert ref==-Fraction(1,414720)
    assert evaluate(folder,executable,graphs[10],source,[(edge,0)])==ref
    print('K4 and wheel current insertions agree with both propagator cuts.',flush=True)


if __name__=='__main__':
    from replay import run
    run(check)
