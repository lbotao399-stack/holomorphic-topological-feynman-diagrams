"""Integrated geometric four-loop kernels in the graph order of four-loop-graphs.json."""
from fractions import Fraction
from functools import cache
from math import prod
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent


@cache
def columns(index):
    if not isinstance(index,int) or not 0<=index<50:
        raise ValueError('the graph index must be an integer from 0 through 49')
    start = next(i for i in (39,11,3,0) if index>=i)
    record = json.loads((HERE/f'kernels/core{start:02}.json').read_text())
    return {tuple(power):tuple(map(Fraction,coeff))
            for power,coeff in zip(record['powers'],record['graphs'][index-start],strict=True)
            if coeff!=['0','0']}


def kernel(index,momenta):
    """Return the multipliers of pi^2 and pi^4 for the four incoming two-component momenta."""
    if len(momenta)!=4 or any(len(p)!=2 for p in momenta) or any(map(sum,zip(*momenta))):
        raise ValueError('four conserved pairs of holomorphic momenta are required')
    minors=[momenta[i][0]*momenta[j][1]-momenta[i][1]*momenta[j][0]
            for i,j in ((0,1),(0,2),(1,2))]
    values=[0,0]
    for powers,coefficients in columns(index).items():
        m=prod(x**j for x,j in zip(minors,powers))
        for j,c in enumerate(coefficients):
            values[j]+=c*m
    return tuple(values)
