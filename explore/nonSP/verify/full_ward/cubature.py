"""Dual time interpolation for the four-loop graph period density."""
from fractions import Fraction as F
from functools import cache
from itertools import product
from math import factorial,lcm,prod
from pathlib import Path
import sys

import networkx as nx
import numpy as np
import sympy as sp

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'prism_forcer'))
from catalogue import graph_catalogue,wheel_table
from four_loop import k4orbit,checked_cut,reduced_powers,reduction
from wheel import orbit as wheel_orbit
from kernel import tree_flows
from check import grid_input
from vk import G


def write_input(index,output):
    edges=graph_catalogue()[index]
    graph=G(6,edges)
    source=[[int(i==j) for j in range(5)] for i in range(5)]+[[-1]*5]
    trees=tree_flows(graph,source)
    _,records=grid_input([(0,0)]*6,edges)
    minors=[rec for rec in records if len(rec)==9]
    moves,name,keep=reduction(edges)
    assert name in ('edge','kite','wheel')
    with output.open('w') as out:
        print(len(trees),len(minors),file=out)
        for miss,flows in trees.items():
            print(*miss,*(x for row in flows[1:] for x in row),file=out)
        for row in minors:
            print(*row,file=out)
        values=rule(graph,moves,name,keep)
        for denom,column in values:
            print(denom,*column,file=out)


def rule(graph,moves,name,keep):
    # Strict convergence on each scaling subset. On the other monomials use
    # the zero extension of the functional; the Kirchhoff numerator has no
    # such monomials for any sources.
    constraints=[]
    for mask in range(1,511):
        es=tuple(i for i in range(9) if mask>>i&1)
        g=nx.Graph(graph.edges[e] for e in es)
        h=len(es)-len(g)+nx.number_connected_components(g)
        if len(es)<=6*h:constraints.append((es,6*h-len(es)))
    table=np.zeros((3,)+(4,)*8,dtype=object)
    for rem in product(range(4),repeat=8):
        k=15-sum(rem)
        if not 0<=k<=3:continue
        ex=(k,)+rem
        if not all(sum(ex[e] for e in es)>threshold for es,threshold in constraints):
            continue
        powers,w=reduced_powers(ex,6,moves,keep)
        if name=='edge':
            assert powers==(0,)
            values=(w,)
        else:
            w*=F(prod(factorial(a-1) for a in powers),factorial(5))
            if name=='kite':
                base=checked_cut(k4orbit(powers))
            else:
                base=wheel_table()[wheel_orbit(powers)]
            values=[w*c for c in base]
        for j,v in enumerate(values):
            table[(j,)+rem]=v
    inv=sp.Matrix([[x**k for k in range(4)] for x in range(1,5)]).inv()
    mat=np.array((6*inv.T).tolist(),dtype=np.int64).astype(object)
    rules=[]
    for channel in table:
        denominator=lcm(*(v.denominator for v in channel.flat))
        arr=np.array([int(v*denominator) for v in channel.flat],dtype=object).reshape(channel.shape)
        for axis in range(8):
            view=np.moveaxis(arr,axis,0)
            arr=np.moveaxis((mat@view.reshape(4,-1)).reshape(view.shape),0,axis)
        rules.append((denominator*6**8,arr.flatten()))
    return rules
